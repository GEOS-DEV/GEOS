#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-2.1-only

"""Fetch an exact-base artifact published by a successful coverage job."""

import argparse
import http.client
import io
import json
import os
import pathlib
import re
import stat
import sys
import tempfile
import urllib.error
import urllib.parse
import urllib.request
import zipfile


SCRIPT_DIRECTORY = pathlib.Path( __file__ ).resolve().parent
if str( SCRIPT_DIRECTORY ) not in sys.path:
    sys.path.insert( 0, str( SCRIPT_DIRECTORY ) )
import check_coverage_thresholds as coverage_contract
import zlib


API_ROOT = "https://api.github.com"
API_VERSION = "2022-11-28"
USER_AGENT = "geos-coverage-baseline-fetch/1"
DOWNLOAD_LIMIT = 10 * 1024 * 1024
JSON_RESPONSE_LIMIT = 2 * 1024 * 1024
PER_PAGE = 100
MAX_PAGES = 10
REQUEST_TIMEOUT_SECONDS = 30
READ_CHUNK_SIZE = 64 * 1024
EXPECTED_FILES = frozenset(
    ( "coverage-summary.json", "llvm-summary.json", "toolchain.txt" )
)
REQUIRED_FILES = frozenset( ( "coverage-summary.json", ) )
STATUS_FILE = "coverage-baseline-status.json"
SHA_RE = re.compile( r"^[0-9a-fA-F]{40}$" )
REPOSITORY_PART_RE = re.compile( r"^[A-Za-z0-9_.-]+$" )
WORKFLOW_RE = re.compile( r"^[A-Za-z0-9_.-]+$" )
ENVIRONMENT_NAME_RE = re.compile( r"^[A-Za-z_][A-Za-z0-9_]*$" )
TOKEN_RE = re.compile( r"^[\x21-\x7e]+$" )
ALLOWED_ZIP_CONTENT_TYPES = frozenset(
    ( "application/octet-stream", "application/x-zip-compressed", "application/zip" )
)
SENSITIVE_REDIRECT_HEADERS = frozenset(
    (
        "authorization",
        "cookie",
        "cookie2",
        "proxy-authorization",
        "x-github-api-version",
    )
)


class BaselineFetchError( Exception ):
    """A malformed input, unsafe artifact, or GitHub API failure."""


class BaselineMissing( Exception ):
    """The exact-base run or artifact is not available."""

    def __init__( self, reason_code: str, message: str ) -> None:
        super().__init__( message )
        self.reason_code = reason_code


def _https_url_parts( url: str, description: str ):
    if not isinstance( url, str ) or any(
        ord( character ) < 32 or ord( character ) == 127 for character in url
    ):
        raise BaselineFetchError( f"{description} is not a safe HTTPS URL" )
    try:
        parsed = urllib.parse.urlsplit( url )
        hostname = parsed.hostname
        port = parsed.port
        username = parsed.username
        password = parsed.password
    except ( TypeError, ValueError ) as error:
        raise BaselineFetchError( f"{description} is not a safe HTTPS URL" ) from error
    if (
        parsed.scheme != "https"
        or not hostname
        or username is not None
        or password is not None
    ):
        raise BaselineFetchError( f"{description} is not a safe HTTPS URL" )
    return parsed, hostname.lower(), 443 if port is None else port


class SafeRedirectHandler( urllib.request.HTTPRedirectHandler ):
    """Allow HTTPS redirects without forwarding credentials across origins."""

    def redirect_request( self, request, file_pointer, code, message, headers, url ):
        del file_pointer, code, message, headers
        try:
            redirected_url = urllib.parse.urljoin( request.full_url, url )
        except ( TypeError, ValueError ) as error:
            raise BaselineFetchError(
                "GitHub redirect URL is not a safe HTTPS URL"
            ) from error
        old, old_hostname, old_port = _https_url_parts(
            request.full_url, "GitHub request URL"
        )
        _, new_hostname, new_port = _https_url_parts(
            redirected_url, "GitHub redirect URL"
        )

        redirected_headers = {}
        same_origin = ( old.scheme, old_hostname, old_port ) == (
            "https",
            new_hostname,
            new_port,
        )
        for name, value in request.header_items():
            if not same_origin and name.lower() in SENSITIVE_REDIRECT_HEADERS:
                continue
            redirected_headers[name] = value
        return urllib.request.Request(
            redirected_url,
            headers=redirected_headers,
            origin_req_host=request.origin_req_host,
            unverifiable=True,
            method="GET",
        )


def _default_opener():
    return urllib.request.build_opener( SafeRedirectHandler() )


def _require_content_length( headers, limit: int ) -> None:
    value = headers.get( "Content-Length" )
    if value is None:
        return
    if not isinstance( value, str ) or not re.fullmatch( r"[0-9]+", value ):
        raise BaselineFetchError( "HTTP response has an invalid Content-Length" )
    length = int( value, 10 )
    if length > limit:
        raise BaselineFetchError(
            f"HTTP response exceeds the {limit}-byte download limit"
        )


def _read_bounded( response, limit: int ) -> bytes:
    _require_content_length( response.headers, limit )
    chunks = []
    size = 0
    while True:
        chunk = response.read( min( READ_CHUNK_SIZE, limit + 1 - size ) )
        if not isinstance( chunk, bytes ):
            raise BaselineFetchError( "HTTP response returned non-byte content" )
        if not chunk:
            return b"".join( chunks )
        chunks.append( chunk )
        size += len( chunk )
        if size > limit:
            raise BaselineFetchError(
                f"HTTP response exceeds the {limit}-byte download limit"
            )


def _content_type( headers ) -> str:
    value = headers.get( "Content-Type", "" )
    return value.split( ";", maxsplit=1 )[0].strip().lower()


class GitHubClient:
    def __init__( self, token: str, opener=None ) -> None:
        if not isinstance( token, str ) or not TOKEN_RE.fullmatch( token ):
            raise BaselineFetchError( "GitHub token is empty or invalid" )
        self._token = token
        self._opener = opener if opener is not None else _default_opener()

    def _request( self, path: str, query=None ):
        if not path.startswith( "/" ):
            raise BaselineFetchError( "GitHub API path must be absolute" )
        url = f"{API_ROOT}{path}"
        if query:
            url = f"{url}?{urllib.parse.urlencode( query )}"
        request = urllib.request.Request(
            url,
            headers={
                "Accept": "application/vnd.github+json",
                "Authorization": f"Bearer {self._token}",
                "User-Agent": USER_AGENT,
                "X-GitHub-Api-Version": API_VERSION,
            },
            method="GET",
        )
        try:
            return self._opener.open( request, timeout=REQUEST_TIMEOUT_SECONDS )
        except urllib.error.HTTPError:
            raise
        except BaselineFetchError:
            raise
        except (
            OSError,
            urllib.error.URLError,
            http.client.HTTPException,
        ) as error:
            raise BaselineFetchError( f"GitHub request failed: {error}" ) from error

    def get_json( self, path: str, query=None ) -> dict:
        try:
            response = self._request( path, query )
            with response:
                status = getattr( response, "status", response.getcode() )
                if status != 200:
                    raise BaselineFetchError(
                        f"GitHub JSON request returned HTTP {status}"
                    )
                _, final_hostname, final_port = _https_url_parts(
                    response.geturl(), "GitHub JSON response URL"
                )
                if final_hostname != "api.github.com" or final_port != 443:
                    raise BaselineFetchError( "GitHub JSON request left api.github.com" )
                if _content_type( response.headers ) not in (
                    "application/json",
                    "application/vnd.github+json",
                ):
                    raise BaselineFetchError(
                        "GitHub JSON response has an unexpected Content-Type"
                    )
                payload = _read_bounded( response, JSON_RESPONSE_LIMIT )
        except urllib.error.HTTPError as error:
            raise BaselineFetchError(
                f"GitHub JSON request returned HTTP {error.code}"
            ) from error
        except BaselineFetchError:
            raise
        except (
            OSError,
            urllib.error.URLError,
            http.client.HTTPException,
        ) as error:
            raise BaselineFetchError(
                f"GitHub JSON response failed: {error}"
            ) from error

        try:
            document = json.loads( payload.decode( "utf-8" ) )
        except ( UnicodeDecodeError, ValueError, RecursionError ) as error:
            raise BaselineFetchError( "GitHub returned malformed JSON" ) from error
        if not isinstance( document, dict ):
            raise BaselineFetchError( "GitHub JSON response must be an object" )
        return document

    def download_artifact( self, repository_path: str, artifact_id: int ) -> bytes:
        path = f"{repository_path}/actions/artifacts/{artifact_id}/zip"
        try:
            response = self._request( path )
            with response:
                status = getattr( response, "status", response.getcode() )
                if status in ( 404, 410 ):
                    raise BaselineMissing(
                        "artifact_download_missing",
                        "the exact-base coverage artifact disappeared or expired",
                    )
                if status != 200:
                    raise BaselineFetchError(
                        f"artifact download returned HTTP {status}"
                    )
                _https_url_parts( response.geturl(), "artifact download URL" )
                if _content_type( response.headers ) not in ALLOWED_ZIP_CONTENT_TYPES:
                    raise BaselineFetchError(
                        "artifact download has an unexpected Content-Type"
                    )
                return _read_bounded( response, DOWNLOAD_LIMIT )
        except urllib.error.HTTPError as error:
            if error.code in ( 404, 410 ):
                raise BaselineMissing(
                    "artifact_download_missing",
                    "the exact-base coverage artifact disappeared or expired",
                ) from error
            raise BaselineFetchError(
                f"artifact download returned HTTP {error.code}"
            ) from error
        except BaselineFetchError:
            raise
        except (
            OSError,
            urllib.error.URLError,
            http.client.HTTPException,
        ) as error:
            raise BaselineFetchError(
                f"artifact download response failed: {error}"
            ) from error


def validate_repository( repository: str ) -> tuple[str, str]:
    parts = repository.split( "/" )
    if (
        len( parts ) != 2
        or not parts[0]
        or not parts[1]
        or any( part in ( ".", ".." ) for part in parts )
        or not all( REPOSITORY_PART_RE.fullmatch( part ) for part in parts )
    ):
        raise BaselineFetchError( "repository must have the form owner/repo" )
    return parts[0], parts[1]


def validate_workflow( workflow: str ) -> str:
    if workflow in ( ".", ".." ) or not WORKFLOW_RE.fullmatch( workflow ):
        raise BaselineFetchError( "workflow must be a workflow filename or numeric ID" )
    return workflow


def validate_sha( revision: str, description: str = "base SHA" ) -> str:
    if not isinstance( revision, str ) or not SHA_RE.fullmatch( revision ):
        raise BaselineFetchError( f"{description} must be exactly 40 hexadecimal digits" )
    return revision.lower()


def validate_artifact_prefix( prefix: str ) -> str:
    if (
        not prefix
        or len( prefix ) > 160
        or any( not character.isprintable() for character in prefix )
        or "/" in prefix
        or "\\" in prefix
    ):
        raise BaselineFetchError( "artifact name prefix is invalid" )
    return prefix


def repository_api_path( repository: str ) -> str:
    owner, name = validate_repository( repository )
    return "/repos/{}/{}".format(
        urllib.parse.quote( owner, safe="" ), urllib.parse.quote( name, safe="" )
    )


def paginated_items(
    client: GitHubClient,
    path: str,
    list_key: str,
    query: dict,
):
    expected_total = None
    yielded = 0
    for page in range( 1, MAX_PAGES + 1 ):
        page_query = dict( query )
        page_query.update( { "page": page, "per_page": PER_PAGE } )
        document = client.get_json( path, page_query )
        total = document.get( "total_count" )
        if isinstance( total, bool ) or not isinstance( total, int ) or total < 0:
            raise BaselineFetchError( "GitHub paginated response has invalid total_count" )
        if expected_total is None:
            expected_total = total
        elif total != expected_total:
            raise BaselineFetchError( "GitHub paginated total_count changed between pages" )
        items = document.get( list_key )
        if not isinstance( items, list ):
            raise BaselineFetchError(
                f"GitHub paginated response is missing the {list_key} array"
            )
        if len( items ) > PER_PAGE:
            raise BaselineFetchError(
                "GitHub paginated response exceeds the requested page size"
            )
        yielded += len( items )
        if yielded > expected_total:
            raise BaselineFetchError(
                "GitHub paginated response contains more entries than total_count"
            )
        for item in items:
            yield item
        if yielded == expected_total:
            return
        if len( items ) < PER_PAGE:
            raise BaselineFetchError(
                "GitHub paginated response ended before total_count"
            )
    raise BaselineFetchError( "GitHub paginated response exceeds the page limit" )


def find_exact_run(
    client: GitHubClient, repository_path: str, workflow: str, base_sha: str
) -> dict:
    workflow_component = urllib.parse.quote( validate_workflow( workflow ), safe="" )
    path = f"{repository_path}/actions/workflows/{workflow_component}/runs"
    query = {
        "branch": "develop",
        "event": "push",
        "head_sha": base_sha,
        # The exact baseline artifact is uploaded only after the coverage job
        # succeeds. The monolithic workflow may still fail for an unrelated
        # CUDA/integrated job, so workflow-level conclusion is not the trust
        # boundary for this artifact.
        "status": "completed",
    }
    matches = []
    for run in paginated_items( client, path, "workflow_runs", query ):
        if not isinstance( run, dict ):
            raise BaselineFetchError( "GitHub workflow run entry must be an object" )
        head_sha = run.get( "head_sha" )
        if not isinstance( head_sha, str ) or not SHA_RE.fullmatch( head_sha ):
            raise BaselineFetchError( "GitHub workflow run has an invalid head_sha" )
        if head_sha.lower() != base_sha:
            continue
        if (
            run.get( "event" ) != "push"
            or run.get( "head_branch" ) != "develop"
            or run.get( "status" ) != "completed"
        ):
            continue
        run_id = run.get( "id" )
        attempt = run.get( "run_attempt", 1 )
        if (
            isinstance( run_id, bool )
            or not isinstance( run_id, int )
            or run_id < 1
            or isinstance( attempt, bool )
            or not isinstance( attempt, int )
            or attempt < 1
        ):
            raise BaselineFetchError( "GitHub workflow run has an invalid identity" )
        matches.append( ( attempt, run_id, run ) )

    if not matches:
        raise BaselineMissing(
            "workflow_run_missing",
            f"no completed develop push run exists for exact base {base_sha}",
        )
    matches.sort( key=lambda candidate: ( candidate[1], candidate[0] ), reverse=True )
    return matches[0][2]


def find_exact_artifact(
    client: GitHubClient,
    repository_path: str,
    run_id: int,
    expected_name: str,
) -> dict:
    path = f"{repository_path}/actions/runs/{run_id}/artifacts"
    exact = []
    for artifact in paginated_items( client, path, "artifacts", {} ):
        if not isinstance( artifact, dict ):
            raise BaselineFetchError( "GitHub artifact entry must be an object" )
        if artifact.get( "name" ) != expected_name:
            continue
        expired = artifact.get( "expired" )
        if not isinstance( expired, bool ):
            raise BaselineFetchError( "GitHub artifact has an invalid expired value" )
        artifact_id = artifact.get( "id" )
        if (
            isinstance( artifact_id, bool )
            or not isinstance( artifact_id, int )
            or artifact_id < 1
        ):
            raise BaselineFetchError( "GitHub artifact has an invalid ID" )
        exact.append( artifact )

    available = [artifact for artifact in exact if not artifact["expired"]]
    if len( available ) > 1:
        raise BaselineFetchError( "multiple unexpired exact-base artifacts were found" )
    if available:
        return available[0]
    if exact:
        raise BaselineMissing(
            "artifact_expired", f"exact-base artifact {expected_name} has expired"
        )
    raise BaselineMissing(
        "artifact_missing", f"exact-base artifact {expected_name} was not found"
    )


def _is_regular_zip_entry( entry: zipfile.ZipInfo ) -> bool:
    if entry.is_dir():
        return False
    if entry.create_system == 3:
        mode = entry.external_attr >> 16
        file_type = stat.S_IFMT( mode )
        return file_type in ( 0, stat.S_IFREG )
    if entry.create_system == 0:
        dos_attributes = entry.external_attr & 0xFF
        return not dos_attributes & 0x18  # Neither a volume label nor a directory.
    return False


def validate_summary_revision( payload: bytes, base_sha: str ) -> None:
    try:
        document = json.loads( payload.decode( "utf-8" ) )
    except ( UnicodeDecodeError, ValueError, RecursionError ) as error:
        raise BaselineFetchError( "coverage-summary.json is not valid UTF-8 JSON" ) from error
    if not isinstance( document, dict ):
        raise BaselineFetchError( "coverage-summary.json must contain an object" )
    schema_version = document.get( "schema_version" )
    if (
        isinstance( schema_version, bool )
        or not isinstance( schema_version, int )
        or schema_version < 1
    ):
        raise BaselineFetchError( "coverage-summary.json has an invalid schema_version" )

    if schema_version == coverage_contract.COVERAGE_SUMMARY_SCHEMA_VERSION:
        try:
            document = coverage_contract.validate_summary( document )
        except ValueError as error:
            raise BaselineFetchError(
                f"coverage-summary.json schema 3 contract is invalid: {error}"
            ) from error
    elif schema_version > coverage_contract.COVERAGE_SUMMARY_SCHEMA_VERSION:
        raise BaselineFetchError(
            f"coverage-summary.json schema {schema_version} is unsupported"
        )

    revisions = []
    measurement = document.get( "measurement" )
    if schema_version >= 3:
        if not isinstance( measurement, dict ) or "commit_sha" not in measurement:
            raise BaselineFetchError(
                "coverage-summary.json schema 3+ requires measurement.commit_sha"
            )
    if measurement is not None:
        if not isinstance( measurement, dict ):
            raise BaselineFetchError(
                "coverage-summary.json measurement must be an object"
            )
        if "commit_sha" in measurement:
            revisions.append( ( "measurement.commit_sha", measurement["commit_sha"] ) )

    source = document.get( "source" )
    if source is not None:
        if not isinstance( source, dict ):
            raise BaselineFetchError( "coverage-summary.json source must be an object" )
        for key in ( "commit", "revision" ):
            if key in source:
                revisions.append( ( f"source.{key}", source[key] ) )
    for key in ( "source_revision", "revision" ):
        if key in document:
            revisions.append( ( key, document[key] ) )

    for description, revision in revisions:
        validated = validate_sha( revision, f"coverage summary {description}" )
        if validated != base_sha:
            raise BaselineFetchError(
                f"coverage summary {description} does not match the requested base SHA"
            )


def validate_optional_json( filename: str, payload: bytes ) -> None:
    try:
        document = json.loads( payload.decode( "utf-8" ) )
    except ( UnicodeDecodeError, ValueError, RecursionError ) as error:
        raise BaselineFetchError( f"{filename} is not valid UTF-8 JSON" ) from error
    if not isinstance( document, dict ):
        raise BaselineFetchError( f"{filename} must contain a JSON object" )


def read_safe_archive( archive: bytes, base_sha: str ) -> dict[str, bytes]:
    if len( archive ) > DOWNLOAD_LIMIT:
        raise BaselineFetchError( "coverage artifact archive exceeds the download limit" )
    try:
        with zipfile.ZipFile( io.BytesIO( archive ) ) as bundle:
            entries = bundle.infolist()
            if not entries:
                raise BaselineFetchError( "coverage artifact ZIP is empty" )
            if len( entries ) > len( EXPECTED_FILES ):
                raise BaselineFetchError(
                    "coverage artifact ZIP contains too many entries"
                )
            result = {}
            total_size = 0
            for entry in entries:
                name = entry.filename
                if (
                    not isinstance( name, str )
                    or entry.orig_filename != name
                    or name not in EXPECTED_FILES
                    or pathlib.PurePosixPath( name ).parts != ( name, )
                    or "\\" in name
                ):
                    raise BaselineFetchError(
                        f"coverage artifact contains unexpected entry {name!r}"
                    )
                if name in result:
                    raise BaselineFetchError(
                        f"coverage artifact contains duplicate entry {name}"
                    )
                if not _is_regular_zip_entry( entry ):
                    raise BaselineFetchError(
                        f"coverage artifact entry {name} is not a regular file"
                    )
                if entry.flag_bits & 0x1:
                    raise BaselineFetchError(
                        f"coverage artifact entry {name} is encrypted"
                    )
                if entry.compress_type not in ( zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED ):
                    raise BaselineFetchError(
                        f"coverage artifact entry {name} uses an unsupported compression"
                    )
                if entry.file_size < 0 or entry.file_size > DOWNLOAD_LIMIT:
                    raise BaselineFetchError(
                        f"coverage artifact entry {name} exceeds the extraction limit"
                    )
                total_size += entry.file_size
                if total_size > DOWNLOAD_LIMIT:
                    raise BaselineFetchError(
                        "coverage artifact uncompressed content exceeds the extraction limit"
                    )
                with bundle.open( entry, "r" ) as stream:
                    payload = stream.read( entry.file_size + 1 )
                if len( payload ) != entry.file_size:
                    raise BaselineFetchError(
                        f"coverage artifact entry {name} has an inconsistent size"
                    )
                result[name] = payload
    except BaselineFetchError:
        raise
    except (
        OSError,
        EOFError,
        RuntimeError,
        ValueError,
        zipfile.BadZipFile,
        zipfile.LargeZipFile,
        zlib.error,
    ) as error:
        raise BaselineFetchError( f"coverage artifact is not a safe ZIP: {error}" ) from error

    missing = REQUIRED_FILES - result.keys()
    if missing:
        raise BaselineFetchError(
            "coverage artifact is missing required file coverage-summary.json"
        )
    validate_summary_revision( result["coverage-summary.json"], base_sha )
    if "llvm-summary.json" in result:
        validate_optional_json( "llvm-summary.json", result["llvm-summary.json"] )
    return result


def prepare_output_directory( output_dir: pathlib.Path ) -> pathlib.Path:
    try:
        if output_dir.is_symlink():
            raise BaselineFetchError( "output directory cannot be a symbolic link" )
        output_dir.mkdir( parents=True, exist_ok=True )
        if output_dir.is_symlink() or not output_dir.is_dir():
            raise BaselineFetchError( "output directory is not a safe directory" )
        output_dir = output_dir.resolve( strict=True )
        for filename in ( *EXPECTED_FILES, STATUS_FILE ):
            path = output_dir / filename
            if os.path.lexists( path ):
                if path.is_dir() and not path.is_symlink():
                    raise BaselineFetchError(
                        f"cannot replace directory at output path {filename}"
                    )
                path.unlink()
    except BaselineFetchError:
        raise
    except OSError as error:
        raise BaselineFetchError( f"cannot prepare output directory: {error}" ) from error
    return output_dir


def atomic_write( path: pathlib.Path, payload: bytes ) -> None:
    temporary_path = None
    try:
        descriptor, temporary_name = tempfile.mkstemp(
            prefix=f".{path.name}.", dir=path.parent
        )
        temporary_path = pathlib.Path( temporary_name )
        with os.fdopen( descriptor, "wb" ) as stream:
            stream.write( payload )
            stream.flush()
            os.fsync( stream.fileno() )
        os.chmod( temporary_path, 0o644 )
        os.replace( temporary_path, path )
        temporary_path = None
    except OSError as error:
        raise BaselineFetchError( f"cannot write {path}: {error}" ) from error
    finally:
        if temporary_path is not None:
            try:
                temporary_path.unlink( missing_ok=True )
            except OSError:
                pass


def write_status( output_dir: pathlib.Path, document: dict ) -> None:
    payload = ( json.dumps( document, indent=2, sort_keys=True ) + "\n" ).encode(
        "utf-8"
    )
    atomic_write( output_dir / STATUS_FILE, payload )


def fetch_baseline(
    client: GitHubClient,
    repository: str,
    workflow: str,
    base_sha: str,
    artifact_prefix: str,
) -> tuple[dict[str, bytes], dict]:
    repository_path = repository_api_path( repository )
    workflow = validate_workflow( workflow )
    base_sha = validate_sha( base_sha )
    artifact_prefix = validate_artifact_prefix( artifact_prefix )
    artifact_name = f"{artifact_prefix}{base_sha}"
    if len( artifact_name ) > 255:
        raise BaselineFetchError( "exact artifact name is too long" )

    run = find_exact_run( client, repository_path, workflow, base_sha )
    artifact = find_exact_artifact(
        client, repository_path, run["id"], artifact_name
    )
    archive = client.download_artifact( repository_path, artifact["id"] )
    files = read_safe_archive( archive, base_sha )
    metadata = {
        "artifact_id": artifact["id"],
        "artifact_name": artifact_name,
        "base_sha": base_sha,
        "repository": repository,
        "run_id": run["id"],
        "schema_version": 1,
        "status": "found",
        "workflow": workflow,
    }
    return files, metadata


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Fetch the exact develop-base LLVM coverage artifact"
    )
    parser.add_argument( "--repository", required=True )
    parser.add_argument( "--workflow", required=True )
    parser.add_argument( "--base-sha", required=True )
    parser.add_argument( "--artifact-name-prefix", required=True )
    parser.add_argument( "--token-env", default="GITHUB_TOKEN" )
    parser.add_argument( "--output-dir", required=True, type=pathlib.Path )
    return parser


def main( argv=None, environ=None, opener=None ) -> int:
    args = build_parser().parse_args( argv )
    environment = os.environ if environ is None else environ
    output_dir = None
    try:
        validate_repository( args.repository )
        validate_workflow( args.workflow )
        base_sha = validate_sha( args.base_sha )
        artifact_prefix = validate_artifact_prefix( args.artifact_name_prefix )
        if not ENVIRONMENT_NAME_RE.fullmatch( args.token_env ):
            raise BaselineFetchError( "token environment variable name is invalid" )
        token = environment.get( args.token_env, "" )
        if not isinstance( token, str ) or not TOKEN_RE.fullmatch( token ):
            raise BaselineFetchError(
                f"GitHub token environment variable {args.token_env} "
                "is missing or empty, or contains unsafe characters"
            )
        output_dir = prepare_output_directory( args.output_dir )
        client = GitHubClient( token, opener=opener )
        files, status = fetch_baseline(
            client,
            args.repository,
            args.workflow,
            base_sha,
            artifact_prefix,
        )
        for filename in sorted( files ):
            atomic_write( output_dir / filename, files[filename] )
        write_status( output_dir, status )
        print(
            f"Fetched {status['artifact_name']} from workflow run {status['run_id']}."
        )
        return 0
    except BaselineMissing as error:
        if output_dir is not None:
            try:
                write_status(
                    output_dir,
                    {
                        "artifact_name": f"{args.artifact_name_prefix}{base_sha}",
                        "base_sha": base_sha,
                        "reason": str( error ),
                        "reason_code": error.reason_code,
                        "repository": args.repository,
                        "schema_version": 1,
                        "status": "missing",
                        "workflow": args.workflow,
                    },
                )
            except BaselineFetchError as status_error:
                print( f"coverage baseline error: {status_error}", file=sys.stderr )
                return 2
        print( f"coverage baseline unavailable: {error}", file=sys.stderr )
        return 3
    except BaselineFetchError as error:
        print( f"coverage baseline error: {error}", file=sys.stderr )
        return 2


if __name__ == "__main__":
    sys.exit( main() )
