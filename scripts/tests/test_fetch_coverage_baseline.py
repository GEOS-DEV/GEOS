#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-2.1-only

import contextlib
import importlib.util
import io
import json
import pathlib
import stat
import tempfile
import unittest
import urllib.error
import urllib.parse
import urllib.request
import warnings
import zipfile
from unittest import mock


SCRIPT = pathlib.Path( __file__ ).parents[1] / "fetch_coverage_baseline.py"
SPEC = importlib.util.spec_from_file_location( "coverage_baseline", SCRIPT )
coverage_baseline = importlib.util.module_from_spec( SPEC )
SPEC.loader.exec_module( coverage_baseline )
CHECKER_TEST_SCRIPT = pathlib.Path( __file__ ).with_name(
    "test_check_coverage_thresholds.py"
)
CHECKER_TEST_SPEC = importlib.util.spec_from_file_location(
    "coverage_threshold_fixtures", CHECKER_TEST_SCRIPT
)
coverage_threshold_fixtures = importlib.util.module_from_spec( CHECKER_TEST_SPEC )
CHECKER_TEST_SPEC.loader.exec_module( coverage_threshold_fixtures )

BASE_SHA = "0123456789abcdef0123456789abcdef01234567"
OTHER_SHA = "fedcba9876543210fedcba9876543210fedcba98"
PREFIX = "llvm-source-coverage-baseline-"


class FakeResponse:
    def __init__(
        self,
        payload,
        *,
        status=200,
        content_type="application/json",
        content_length=True,
        final_url=None,
        max_read_size=None,
        read_error=None,
    ):
        self._stream = io.BytesIO( payload )
        self.status = status
        self.headers = { "Content-Type": content_type }
        if content_length:
            self.headers["Content-Length"] = str( len( payload ) )
        self.final_url = final_url
        self.max_read_size = max_read_size
        self.read_error = read_error

    def __enter__( self ):
        return self

    def __exit__( self, exception_type, exception, traceback ):
        return False

    def getcode( self ):
        return self.status

    def geturl( self ):
        return self.final_url

    def read( self, size=-1 ):
        if self.read_error is not None:
            raise self.read_error
        if self.max_read_size is not None and (
            size < 0 or size > self.max_read_size
        ):
            size = self.max_read_size
        return self._stream.read( size )


class QueueOpener:
    def __init__( self, responses ):
        self.responses = list( responses )
        self.requests = []

    def open( self, request, timeout ):
        self.requests.append( ( request, timeout ) )
        if not self.responses:
            raise AssertionError( f"unexpected request: {request.full_url}" )
        response = self.responses.pop( 0 )
        if isinstance( response, BaseException ):
            raise response
        if callable( response ):
            response = response( request )
        if response.final_url is None:
            response.final_url = request.full_url
        return response


def json_response( document, **kwargs ):
    return FakeResponse( json.dumps( document ).encode( "utf-8" ), **kwargs )


def workflow_run( run_id=101, sha=BASE_SHA, **overrides ):
    document = {
        "conclusion": "success",
        "event": "push",
        "head_branch": "develop",
        "head_sha": sha,
        "id": run_id,
        "run_attempt": 1,
        "status": "completed",
    }
    document.update( overrides )
    return document


def artifact( artifact_id=202, name=None, expired=False ):
    return {
        "expired": expired,
        "id": artifact_id,
        "name": name if name is not None else f"{PREFIX}{BASE_SHA}",
    }


def coverage_summary( revision=BASE_SHA, include_revision=True ):
    document = (
        coverage_threshold_fixtures.summary()
        if include_revision
        else { "schema_version": 2 }
    )
    if include_revision:
        document["measurement"]["commit_sha"] = revision
    return json.dumps( document, sort_keys=True ).encode( "utf-8" )


def make_zip( entries ):
    output = io.BytesIO()
    with warnings.catch_warnings():
        warnings.simplefilter( "ignore", UserWarning )
        with zipfile.ZipFile( output, "w" ) as bundle:
            for entry in entries:
                if len( entry ) == 2:
                    name, payload = entry
                    bundle.writestr( name, payload, compress_type=zipfile.ZIP_DEFLATED )
                else:
                    name, payload, mode = entry
                    information = zipfile.ZipInfo( name )
                    information.create_system = 3
                    information.external_attr = mode << 16
                    information.compress_type = zipfile.ZIP_DEFLATED
                    bundle.writestr( information, payload )
    return output.getvalue()


def successful_responses( archive, *, artifact_document=None ):
    artifacts = (
        [artifact()] if artifact_document is None else artifact_document
    )
    return [
        json_response( { "total_count": 1, "workflow_runs": [workflow_run()] } ),
        json_response( { "artifacts": artifacts, "total_count": len( artifacts ) } ),
        FakeResponse( archive, content_type="application/zip" ),
    ]


def cli_arguments( output_dir ):
    return [
        "--repository",
        "LLNL/GEOS",
        "--workflow",
        "ci_tests.yml",
        "--base-sha",
        BASE_SHA,
        "--artifact-name-prefix",
        PREFIX,
        "--token-env",
        "TEST_TOKEN",
        "--output-dir",
        str( output_dir ),
    ]


def run_main( opener, output_dir, environment=None ):
    standard_output = io.StringIO()
    standard_error = io.StringIO()
    if environment is None:
        environment = { "TEST_TOKEN": "secret-token" }
    with contextlib.redirect_stdout( standard_output ), contextlib.redirect_stderr(
        standard_error
    ):
        status = coverage_baseline.main(
            cli_arguments( output_dir ), environ=environment, opener=opener
        )
    return status, standard_output.getvalue(), standard_error.getvalue()


class CoverageBaselineFetchTests( unittest.TestCase ):
    def test_success_fetches_only_exact_flat_files( self ):
        archive = make_zip(
            [
                ( "coverage-summary.json", coverage_summary() ),
                ( "llvm-summary.json", b'{"data": []}\n' ),
                ( "toolchain.txt", b"LLVM 20.1.2\n" ),
            ]
        )
        wrong_name = artifact(
            200, name=f"{PREFIX}{BASE_SHA}-extra", expired=False
        )
        opener = QueueOpener(
            successful_responses( archive, artifact_document=[wrong_name, artifact()] )
        )
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_dir = pathlib.Path( temporary_directory )
            ( output_dir / "toolchain.txt" ).write_text( "stale", encoding="utf-8" )
            status, standard_output, standard_error = run_main( opener, output_dir )

            self.assertEqual( status, 0 )
            self.assertEqual( standard_error, "" )
            self.assertIn( f"{PREFIX}{BASE_SHA}", standard_output )
            self.assertEqual(
                json.loads(
                    ( output_dir / "coverage-summary.json" ).read_text(
                        encoding="utf-8"
                    )
                )["measurement"]["commit_sha"],
                BASE_SHA,
            )
            self.assertEqual(
                ( output_dir / "toolchain.txt" ).read_text( encoding="utf-8" ),
                "LLVM 20.1.2\n",
            )
            status_document = json.loads(
                ( output_dir / coverage_baseline.STATUS_FILE ).read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual( status_document["status"], "found" )
            self.assertEqual( status_document["run_id"], 101 )
            self.assertEqual( status_document["artifact_id"], 202 )

        self.assertEqual( len( opener.requests ), 3 )
        runs_url = urllib.parse.urlsplit( opener.requests[0][0].full_url )
        query = urllib.parse.parse_qs( runs_url.query )
        self.assertEqual( query["head_sha"], [BASE_SHA] )
        self.assertEqual( query["branch"], ["develop"] )
        self.assertEqual( query["event"], ["push"] )
        self.assertEqual( query["status"], ["completed"] )
        self.assertEqual(
            opener.requests[2][0].full_url,
            "https://api.github.com/repos/LLNL/GEOS/actions/artifacts/202/zip",
        )
        self.assertEqual(
            opener.requests[0][0].get_header( "Authorization" ),
            "Bearer secret-token",
        )

    def test_schema_without_embedded_revision_is_accepted( self ):
        archive = make_zip(
            [
                (
                    "coverage-summary.json",
                    coverage_summary( include_revision=False ),
                )
            ]
        )
        files = coverage_baseline.read_safe_archive( archive, BASE_SHA )
        self.assertEqual( set( files ), { "coverage-summary.json" } )

    def test_schema_three_requires_embedded_revision( self ):
        document = coverage_threshold_fixtures.summary()
        del document["measurement"]
        archive = make_zip(
            [( "coverage-summary.json", json.dumps( document ).encode( "utf-8" ) )]
        )
        with self.assertRaisesRegex(
            coverage_baseline.BaselineFetchError, "measurement"
        ):
            coverage_baseline.read_safe_archive( archive, BASE_SHA )

    def test_missing_run_is_exit_three_and_clears_stale_outputs( self ):
        opener = QueueOpener(
            [json_response( { "total_count": 0, "workflow_runs": [] } )]
        )
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_dir = pathlib.Path( temporary_directory )
            ( output_dir / "coverage-summary.json" ).write_text(
                "stale", encoding="utf-8"
            )
            status, _, standard_error = run_main( opener, output_dir )
            self.assertEqual( status, 3 )
            self.assertIn( "unavailable", standard_error )
            self.assertFalse( ( output_dir / "coverage-summary.json" ).exists() )
            reason = json.loads(
                ( output_dir / coverage_baseline.STATUS_FILE ).read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual( reason["status"], "missing" )
            self.assertEqual( reason["reason_code"], "workflow_run_missing" )
            self.assertEqual( reason["base_sha"], BASE_SHA )

    def test_expired_or_inexact_artifact_is_exit_three( self ):
        cases = (
            (
                [artifact( expired=True )],
                "artifact_expired",
            ),
            (
                [artifact( name=f"{PREFIX}{BASE_SHA}-not-exact" )],
                "artifact_missing",
            ),
        )
        for artifacts, reason_code in cases:
            with self.subTest( reason_code=reason_code ):
                opener = QueueOpener(
                    [
                        json_response(
                            { "total_count": 1, "workflow_runs": [workflow_run()] }
                        ),
                        json_response(
                            { "artifacts": artifacts, "total_count": len( artifacts ) }
                        ),
                    ]
                )
                with tempfile.TemporaryDirectory() as temporary_directory:
                    output_dir = pathlib.Path( temporary_directory )
                    status, _, _ = run_main( opener, output_dir )
                    self.assertEqual( status, 3 )
                    reason = json.loads(
                        ( output_dir / coverage_baseline.STATUS_FILE ).read_text(
                            encoding="utf-8"
                        )
                    )
                    self.assertEqual( reason["reason_code"], reason_code )

    def test_download_404_is_a_missing_baseline( self ):
        download_url = (
            "https://api.github.com/repos/LLNL/GEOS/actions/artifacts/202/zip"
        )
        missing_responses = (
            urllib.error.HTTPError( download_url, 404, "missing", {}, None ),
            FakeResponse( b"", status=410, content_type="application/zip" ),
        )
        for missing_response in missing_responses:
            with self.subTest( response=missing_response ):
                opener = QueueOpener(
                    [
                        json_response(
                            {
                                "total_count": 1,
                                "workflow_runs": [workflow_run()],
                            }
                        ),
                        json_response(
                            { "artifacts": [artifact()], "total_count": 1 }
                        ),
                        missing_response,
                    ]
                )
                with tempfile.TemporaryDirectory() as temporary_directory:
                    output_dir = pathlib.Path( temporary_directory )
                    status, _, _ = run_main( opener, output_dir )
                    self.assertEqual( status, 3 )
                    reason = json.loads(
                        ( output_dir / coverage_baseline.STATUS_FILE ).read_text(
                            encoding="utf-8"
                        )
                    )
                    self.assertEqual(
                        reason["reason_code"], "artifact_download_missing"
                    )

    def test_http_and_json_failures_are_exit_two( self ):
        runs_url = (
            "https://api.github.com/repos/LLNL/GEOS/actions/workflows/"
            "ci_tests.yml/runs"
        )
        cases = (
            QueueOpener(
                [urllib.error.HTTPError( runs_url, 403, "forbidden", {}, None )]
            ),
            QueueOpener(
                [FakeResponse( b"not json", content_type="application/json" )]
            ),
            QueueOpener(
                [json_response( { "workflow_runs": [], "total_count": "zero" } )]
            ),
        )
        for opener in cases:
            with self.subTest( opener=opener ):
                with tempfile.TemporaryDirectory() as temporary_directory:
                    status, _, standard_error = run_main(
                        opener, pathlib.Path( temporary_directory )
                    )
                    self.assertEqual( status, 2 )
                    self.assertIn( "coverage baseline error", standard_error )

    def test_response_content_type_final_url_and_size_are_strict( self ):
        document = { "total_count": 0, "workflow_runs": [] }
        cases = (
            json_response( document, content_type="text/plain" ),
            json_response( document, final_url="https://example.com/runs" ),
            json_response(
                document,
                final_url="https://api.github.com:invalid/runs",
            ),
            FakeResponse(
                b"{}",
                content_type="application/json",
                content_length=False,
            ),
        )
        cases[3].headers["Content-Length"] = str(
            coverage_baseline.JSON_RESPONSE_LIMIT + 1
        )
        for response in cases:
            with self.subTest( response=response ):
                client = coverage_baseline.GitHubClient(
                    "token", opener=QueueOpener( [response] )
                )
                with self.assertRaises( coverage_baseline.BaselineFetchError ):
                    client.get_json( "/repos/LLNL/GEOS/actions/runs" )

    def test_partial_reads_are_completed_and_read_failures_are_wrapped( self ):
        document = { "total_count": 0, "workflow_runs": [] }
        client = coverage_baseline.GitHubClient(
            "token",
            opener=QueueOpener(
                [json_response( document, max_read_size=3 )]
            ),
        )
        self.assertEqual(
            client.get_json( "/repos/LLNL/GEOS/actions/runs" ), document
        )

        client = coverage_baseline.GitHubClient(
            "token",
            opener=QueueOpener(
                [FakeResponse( b"", read_error=OSError( "read failed" ) )]
            ),
        )
        with self.assertRaisesRegex(
            coverage_baseline.BaselineFetchError, "response failed"
        ):
            client.get_json( "/repos/LLNL/GEOS/actions/runs" )

    def test_latest_exact_completed_develop_push_is_selected( self ):
        runs = [
            workflow_run( 100, conclusion="failure" ),
            workflow_run( 101, head_branch="other" ),
            workflow_run( 102, sha=OTHER_SHA ),
            workflow_run( 103, run_attempt=2 ),
            workflow_run( 104, run_attempt=1 ),
        ]
        client = coverage_baseline.GitHubClient(
            "token",
            opener=QueueOpener(
                [json_response( { "total_count": len( runs ), "workflow_runs": runs } )]
            ),
        )
        selected = coverage_baseline.find_exact_run(
            client, "/repos/LLNL/GEOS", "ci_tests.yml", BASE_SHA
        )
        self.assertEqual( selected["id"], 104 )

    def test_pagination_finds_an_exact_run_on_a_later_page( self ):
        with mock.patch.object( coverage_baseline, "PER_PAGE", 2 ):
            opener = QueueOpener(
                [
                    json_response(
                        {
                            "total_count": 3,
                            "workflow_runs": [
                                workflow_run( 1, sha=OTHER_SHA ),
                                workflow_run( 2, sha=OTHER_SHA ),
                            ],
                        }
                    ),
                    json_response(
                        {
                            "total_count": 3,
                            "workflow_runs": [workflow_run( 3 )],
                        }
                    ),
                ]
            )
            client = coverage_baseline.GitHubClient( "token", opener=opener )
            selected = coverage_baseline.find_exact_run(
                client, "/repos/LLNL/GEOS", "12345", BASE_SHA
            )
        self.assertEqual( selected["id"], 3 )
        second_query = urllib.parse.parse_qs(
            urllib.parse.urlsplit( opener.requests[1][0].full_url ).query
        )
        self.assertEqual( second_query["page"], ["2"] )

    def test_pagination_rejects_a_short_incomplete_response( self ):
        client = coverage_baseline.GitHubClient(
            "token",
            opener=QueueOpener(
                [
                    json_response(
                        {
                            "total_count": 2,
                            "workflow_runs": [workflow_run()],
                        }
                    )
                ]
            ),
        )
        with self.assertRaisesRegex(
            coverage_baseline.BaselineFetchError, "before total_count"
        ):
            coverage_baseline.find_exact_run(
                client, "/repos/LLNL/GEOS", "ci_tests.yml", BASE_SHA
            )

    def test_duplicate_exact_artifacts_are_rejected( self ):
        client = coverage_baseline.GitHubClient(
            "token",
            opener=QueueOpener(
                [
                    json_response(
                        {
                            "artifacts": [artifact( 1 ), artifact( 2 )],
                            "total_count": 2,
                        }
                    )
                ]
            ),
        )
        with self.assertRaisesRegex(
            coverage_baseline.BaselineFetchError, "multiple"
        ):
            coverage_baseline.find_exact_artifact(
                client,
                "/repos/LLNL/GEOS",
                10,
                f"{PREFIX}{BASE_SHA}",
            )

    def test_zip_rejects_unexpected_traversal_duplicate_and_missing_entries( self ):
        valid_summary = coverage_summary()
        archives = (
            make_zip( [( "../coverage-summary.json", valid_summary )] ),
            make_zip(
                [
                    ( "coverage-summary.json", valid_summary ),
                    ( "unexpected.txt", b"unexpected" ),
                ]
            ),
            make_zip(
                [
                    ( "coverage-summary.json", valid_summary ),
                    ( "coverage-summary.json", valid_summary ),
                ]
            ),
            make_zip( [( "toolchain.txt", b"LLVM" )] ),
        )
        for archive in archives:
            with self.subTest():
                with self.assertRaises( coverage_baseline.BaselineFetchError ):
                    coverage_baseline.read_safe_archive( archive, BASE_SHA )

    def test_zip_rejects_symlinks_and_special_files( self ):
        modes = ( stat.S_IFLNK | 0o777, stat.S_IFIFO | 0o600 )
        for mode in modes:
            archive = make_zip(
                [( "coverage-summary.json", b"target", mode )]
            )
            with self.subTest( mode=mode ):
                with self.assertRaisesRegex(
                    coverage_baseline.BaselineFetchError, "regular file"
                ):
                    coverage_baseline.read_safe_archive( archive, BASE_SHA )

        dos_directory = zipfile.ZipInfo( "coverage-summary.json" )
        dos_directory.create_system = 0
        dos_directory.external_attr = 0x10
        self.assertFalse( coverage_baseline._is_regular_zip_entry( dos_directory ) )

        unknown_system = zipfile.ZipInfo( "coverage-summary.json" )
        unknown_system.create_system = 7
        self.assertFalse( coverage_baseline._is_regular_zip_entry( unknown_system ) )

    def test_zip_rejects_invalid_json_revision_mismatch_and_malformed_zip( self ):
        archives = (
            make_zip( [( "coverage-summary.json", b"not json" )] ),
            make_zip(
                [( "coverage-summary.json", coverage_summary( OTHER_SHA ) )]
            ),
            b"not a zip archive",
        )
        for archive in archives:
            with self.subTest():
                with self.assertRaises( coverage_baseline.BaselineFetchError ):
                    coverage_baseline.read_safe_archive( archive, BASE_SHA )

    def test_zip_bomb_is_rejected_by_uncompressed_limit( self ):
        oversized = b" " * ( coverage_baseline.DOWNLOAD_LIMIT + 1 )
        archive = make_zip( [( "coverage-summary.json", oversized )] )
        self.assertLess( len( archive ), coverage_baseline.DOWNLOAD_LIMIT )
        with self.assertRaisesRegex(
            coverage_baseline.BaselineFetchError, "extraction limit"
        ):
            coverage_baseline.read_safe_archive( archive, BASE_SHA )

    def test_download_body_is_limited_even_without_content_length( self ):
        oversized = b"x" * ( coverage_baseline.DOWNLOAD_LIMIT + 1 )
        response = FakeResponse(
            oversized,
            content_type="application/zip",
            content_length=False,
            final_url="https://objects.example.test/artifact.zip",
        )
        client = coverage_baseline.GitHubClient(
            "token", opener=QueueOpener( [response] )
        )
        with self.assertRaisesRegex(
            coverage_baseline.BaselineFetchError, "download limit"
        ):
            client.download_artifact( "/repos/LLNL/GEOS", 1 )

    def test_download_response_type_and_url_are_strict( self ):
        responses = (
            FakeResponse( b"zip", content_type="text/plain" ),
            FakeResponse(
                b"zip",
                content_type="application/zip",
                final_url="http://objects.example.test/artifact.zip",
            ),
        )
        for response in responses:
            with self.subTest( response=response ):
                client = coverage_baseline.GitHubClient(
                    "token", opener=QueueOpener( [response] )
                )
                with self.assertRaises( coverage_baseline.BaselineFetchError ):
                    client.download_artifact( "/repos/LLNL/GEOS", 1 )

    def test_cross_origin_redirect_strips_authorization( self ):
        handler = coverage_baseline.SafeRedirectHandler()
        request = urllib.request.Request(
            "https://api.github.com/repos/LLNL/GEOS/actions/artifacts/1/zip",
            headers={
                "Authorization": "Bearer secret",
                "Cookie": "session=secret",
                "Proxy-Authorization": "Basic secret",
                "X-GitHub-Api-Version": coverage_baseline.API_VERSION,
                "User-Agent": "test",
            },
        )
        redirected = handler.redirect_request(
            request,
            None,
            302,
            "Found",
            {},
            "https://objects.example.test/artifact.zip",
        )
        headers = {
            name.lower(): value for name, value in redirected.header_items()
        }
        self.assertNotIn( "authorization", headers )
        self.assertNotIn( "cookie", headers )
        self.assertNotIn( "proxy-authorization", headers )
        self.assertNotIn( "x-github-api-version", headers )
        self.assertEqual( headers["user-agent"], "test" )
        unsafe_urls = (
            "http://example.test/artifact.zip",
            "https://@objects.example.test/artifact.zip",
            "https://api.github.com:invalid/artifact.zip",
            "https://[malformed/artifact.zip",
        )
        for unsafe_url in unsafe_urls:
            with self.subTest( unsafe_url=unsafe_url ):
                with self.assertRaises( coverage_baseline.BaselineFetchError ):
                    handler.redirect_request(
                        request, None, 302, "Found", {}, unsafe_url
                    )

    def test_unsafe_inputs_and_output_symlink_are_exit_two( self ):
        invalid_argument_sets = (
            ( "--repository", "../GEOS" ),
            ( "--workflow", "../ci_tests.yml" ),
            ( "--base-sha", "short" ),
            ( "--artifact-name-prefix", "bad/name" ),
            ( "--token-env", "BAD-NAME" ),
        )
        for option, value in invalid_argument_sets:
            with self.subTest( option=option ):
                with tempfile.TemporaryDirectory() as temporary_directory:
                    arguments = cli_arguments( pathlib.Path( temporary_directory ) )
                    arguments[arguments.index( option ) + 1] = value
                    with contextlib.redirect_stderr( io.StringIO() ):
                        status = coverage_baseline.main(
                            arguments,
                            environ={ "TEST_TOKEN": "token", "BAD-NAME": "token" },
                            opener=QueueOpener( [] ),
                        )
                    self.assertEqual( status, 2 )

        with tempfile.TemporaryDirectory() as temporary_directory:
            root = pathlib.Path( temporary_directory )
            target = root / "target"
            target.mkdir()
            output_link = root / "output"
            output_link.symlink_to( target, target_is_directory=True )
            with contextlib.redirect_stderr( io.StringIO() ):
                status = coverage_baseline.main(
                    cli_arguments( output_link ),
                    environ={ "TEST_TOKEN": "token" },
                    opener=QueueOpener( [] ),
                )
            self.assertEqual( status, 2 )

    def test_missing_token_is_exit_two_without_an_http_request( self ):
        environments = (
            {},
            { "TEST_TOKEN": "token\nInjected: value" },
            { "TEST_TOKEN": "token\0suffix" },
        )
        for environment in environments:
            with self.subTest( environment=environment ):
                opener = QueueOpener( [] )
                with tempfile.TemporaryDirectory() as temporary_directory:
                    status, _, standard_error = run_main(
                        opener,
                        pathlib.Path( temporary_directory ),
                        environment=environment,
                    )
                self.assertEqual( status, 2 )
                self.assertIn( "missing or empty", standard_error )
                self.assertEqual( opener.requests, [] )


if __name__ == "__main__":
    unittest.main()
