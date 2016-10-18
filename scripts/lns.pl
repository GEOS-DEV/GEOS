#!/usr/bin/env perl

# lns -- create a symbolic link. Alternative to "ln -s".
# This program works more like "cp", in that the source path name is not
# taken literally.
#    ln -s filename /tmp
# creates a link
#    /tmp/filename -> filename
# whereas we would prefer
#    /tmp/filename -> /home/me/filename
# or wherever the file *really* was.
#
# Usage: lns [-afF] file1 file2
#     or lns [-af] file1 [file2...] dir
#
# Where:
#   -a means absolute - "symlink /usr/bin/argh /usr/local/bin/argh" produces
#      a relative link "/usr/bin/argh -> ../local/bin/argh", but using the
#      -a option will give a real absolute link.
#   -f means forceful - overwrite the target filename if it exists *and* is
#      a link. You can't accidentally overwrite real files like this.
#   -q means quiet - don't complain if we fail to do the job.
#   -v means verbose - say what we're doing.
#   -F means FILE - forces interpretation to be the "file1 file2" syntax,
#      even if file2 is a link to a directory. This option implies -f.

use Cwd;
use POSIX; # for opendir and friends

$usage =
  "usage: lns [flags] srcfile destfile\n".
  "   or: lns [flags] srcfile [srcfile...] destdir\n".
  "where: -a               create symlinks with absolute path names\n".
  "       -f               overwrite existing symlink at target location\n".
  "       -F               like -f, but works even if target is link to dir\n".
  "       -e               tolerate an identical symlink at target location\n".
  "       -r               recursively construct a directory tree which\n".
  "                        mirrors the source, with symlinks to all files\n".
  "       -v               verbosely log activity (repeat for more verbosity)\n".
  "       -q               suppress error messages on failure\n".
  " also: lns --version    report version number\n" .
  "       lns --help       display this help text\n" .
  "       lns --licence    display (MIT) licence text\n";

$licence =
  "lns is copyright 1999,2004 Simon Tatham.\n" .
  "\n" .
  "Permission is hereby granted, free of charge, to any person\n" .
  "obtaining a copy of this software and associated documentation files\n" .
  "(the \"Software\"), to deal in the Software without restriction,\n" .
  "including without limitation the rights to use, copy, modify, merge,\n" .
  "publish, distribute, sublicense, and/or sell copies of the Software,\n" .
  "and to permit persons to whom the Software is furnished to do so,\n" .
  "subject to the following conditions:\n" .
  "\n" .
  "The above copyright notice and this permission notice shall be\n" .
  "included in all copies or substantial portions of the Software.\n" .
  "\n" .
  "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND,\n" .
  "EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF\n" .
  "MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND\n" .
  "NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS\n" .
  "BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN\n" .
  "ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN\n" .
  "CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\n" .
  "SOFTWARE.\n";

$abs=$force=$quiet=$verbose=$recurse=$tolerate=$FILE=0;
while ($_=shift @ARGV) {
  last if /^--$/;
  unshift (@ARGV, $_), last unless /^-(.*)/;
  if ($1 eq "-help") {
    print STDERR $usage;
    exit 0;
  } elsif ($1 eq "-version") {
    print "lns, version 20161010.cd9f6e3\n";
    exit 0;
  } elsif ($1 eq "-licence" or $1 eq "-license") {
    print $licence;
    exit 0;
  } else {
    foreach $opt (split //, $1) {
	if ($opt eq "a") { $abs=1; }
	elsif ($opt eq "f") { $force=1; }
	elsif ($opt eq "e") { $tolerate=1; }
	elsif ($opt eq "q") { $quiet=1; }
	elsif ($opt eq "r") { $recurse=1; }
	elsif ($opt eq "v") { $verbose++; }
	elsif ($opt eq "F") { $force=$FILE=1; }
	else { die "lns: unrecognised option '-$1'\n"; }
    }
  }
}

die $usage if $#ARGV < 1;

die "lns: multiple source files specified with -F option\n"
  if $#ARGV > 1 && $FILE;
die "lns: -q (quiet) and -v (verbose) options both specified\n"
  if $quiet && $verbose;

$target = pop @ARGV;
die "lns: multiple source files specified, $target not a directory\n"
  if $#ARGV > 0 && !-d $target;

$multiple = (-d $target && !$FILE);

$target =~ s/// if $target =~ /\/$/;    # strip trailing slash if present

if ($multiple) {
  foreach $source (@ARGV) {
    # We must path-normalise $source _before_ looking for the final
    # filename component, to deal with the case of `lns . subdir'
    # in which we want the link to be called subdir/<dirname> rather
    # than subdir/. .
    $source = &normalise($source);
    $source =~ /^(.*\/)?([^\/]*)$/;     # find final file name component
    &makelink($source, "$target/$2");   # actually make a link
  }
} else {
  $source = $ARGV[0];                   # only one source file
  &makelink($source, $target);          # make the link
}

sub makelink {
  local ($source, $target) = @_;
  # Calculate the absolute path names of both source and target.
  $source = &normalise($source);
  $target = &normalise($target);

  # If we're in Relative mode (the default), calculate the relative path
  # name we will reference the source by.
  $sourcename = $abs ? $source : &relname($source, $target);

  my $donothing = 0;
  my $recursing = $recurse && -d $source;
  my $ok;
  # If the target exists...
  if (-e $target || readlink $target) {
    if ($recursing && -d $target) {
      # If it's a directory and we're in recursive mode, just do nothing
      # and work around it.
      $donothing = 1;
    } elsif ($force && readlink $target) {
      # If it's a symlink and we're in Force mode, remove it and carry on.
      unlink $target || die "lns: unable to remove link $target\n";
      # Report that if in Verbose mode.
      warn "lns: removing existing target link $target\n" if $verbose;
    } elsif ($tolerate && $sourcename eq readlink $target) {
      # If the symlink already exists and is the same one we would
      # have created anyway, and we're in -e mode, do nothing.
      warn "lns: nothing to do, $target already points to $source\n"
          if $verbose;
      return;
    } else {
      # Otherwise, fail. Report that fact if not in Quiet mode.
      warn "lns: failed to link $source to $target: target exists\n"
        if !$quiet;
      return;
    }
  }

  if ($recursing) {
    # Make the directory.
    if ($donothing) {
      warn "lns: directory $target already exists, no need to create it\n"
          if $verbose;
      $ok = 1;
    } else {
      warn "lns: making directory $target\n"
          if $verbose;
      if (mkdir $target) {
        $ok = 1;
      } else {
        warn "lns: unable to make directory '$target': $!\n";
      }
    }
    # Now recurse into it.
    if ($ok) {
      my $dh = POSIX::opendir($source);
      my @files = POSIX::readdir($dh);
      my $f;
      POSIX::closedir($dh);
      foreach $f (@files) {
        next if $f eq "." or $f eq "..";
        &makelink("$source/$f", "$target/$f");
      }
    }
  } else {
    # Make the link.
    warn "lns: linking $source: $target -> $sourcename\n" if $verbose;
    symlink($sourcename, $target) || die "lns: unable to make link to $target\n";
  }
}

sub normalise {
    # Normalise a path into an absolute one containing no . or ..
    # segments.
    local ($_) = @_;

    warn "lns: path normalisation required on $_\n" if $verbose > 2;

    # Make relative paths absolute.
    $_ = getcwd() . "/" . $_ if !/^\//;

    # Remove "." segments.
    1 while s/^(.*)\/\.(\/.*)?$/$1$2/;

    # Remove redundant slashes.
    s/\/+/\//g;

    # Remove a trailing slash if present.
    s/\/$//;

    # Remove ".." segments. This is the hard bit, because a
    # directory segment that's a _symlink_ doesn't do the obvious
    # thing if followed by "..". But we can't just call realpath,
    # because we do want to preserve symlinks where they _don't_
    # interfere with this sort of work. So the algorithm is:
    #
    #  - Repeatedly search for the rightmost `directory/..'
    # 	 fragment.
    #  - When we find it, one of two cases apply.
    # 	  * If the directory before the .. is not a symlink, we can
    # 	    remove both it and the .. from the string.
    # 	  * If it _is_ a symlink, we substitute it for its link
    # 	    text, and loop round again.
    while (/^(.*)\/((\.|\.\.[^\/]+|\.?[^\/\.][^\/]*)\/\.\.)(\/.*)?$/)
    {
	my ($pre, $dir, $frag, $post) = ($1,$2,$3,$4);
	my $log = "  transforming $_ -> ";
	if (-l "$pre/$frag") {
	    my $linktext = readlink "$pre/$frag";
	    if ($linktext =~ /^\//) { # absolute link
		$_ = $linktext;
	    } else { # relative link
		$_ = "$pre/$linktext";
	    }
	    $_ .= "/.." . $post;
	} else {
	    $_ = $pre . $post;
	}
	$_ = "/" if $_ eq ""; # special case
	s/\/+/\//g; # remove redundant slashes again in case link text had any
	$log .= "$_";
	warn "lns: $log\n" if $verbose > 2;
    }

    # The only place where a ".." fragment might still remain is at
    # the very start of the string, and "/.." is defined to be
    # equivalent to "/".
    1 while s/^\/\.\.(\/(.*))?$/\/$2/;

    warn "lns: path normalisation returned $_\n" if $verbose > 2;

    return $_;
}

sub relname {
  local ($source, $target) = @_;
  local $prefix;

  # Strip the last word off the target (the actual file name) to
  # obtain the target _directory_.
  $target =~ s/\/[^\/]*$//;

  # Our starting prefix is empty. We will add one "../" at a time
  # until we find a match.

  while (1) {

    warn "lns: trying prefix '$prefix': looking for $target as prefix of $source\n" if $verbose > 1;

    # If $target is _precisely_ $source, we are done.
    if ($target eq $source) {
      return "." if $prefix eq "";
      $prefix =~ s/\/$//;
      return $prefix;
    }

    # If $target is a prefix of $source, we are done. (No matter what
    # symlinks may exist on the shared common pathname, if we are
    # linking `a/b/c/foo' to `foo' then a simple relative link will
    # work.)
    if (substr($source, 0, 1 + length $target) eq "$target/") {
      warn "lns: found it\n" if $verbose > 1;
      return $prefix . substr($source, 1 + length $target); # skip the slash
    } elsif ($target eq "/") {
      warn "lns: found it\n" if $verbose > 1;
      return $prefix . substr($source, 1); # special case
    }

    # Otherwise, descend to "..".
    $target = &normalise($target . "/..");

    # Now we have replaced $target with a pathname equivalent to
    # `$target/..'. So add a "../" to $prefix, and try matching
    # again.
    $prefix .= "../";
  }
}
