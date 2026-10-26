#!/usr/bin/env bash
###############################################################################
#                                                                             #
#   deb.sh: the Debian package of PBeFS from a release tarball                #
#                                                                             #
#   Usage:  deb.sh <tarball> binary          libpbefs-dev, unsigned           #
#           deb.sh <tarball> <series> [n]    the signed source package of     #
#                                            X.Y.Z-1~<series><n> for a PPA    #
#                                                                             #
#   Unpacks the tarball, with a copy of it as pbefs_<version>.orig.tar.gz,    #
#   adds packaging/debian and runs dpkg-buildpackage. Nothing is uploaded:    #
#   for a PPA it prints the dput. The output goes to build-deb/ in the        #
#   current directory, next to the copy of the tarball.                       #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
#                                                                             #
###############################################################################

set -euo pipefail

tarball=$(realpath "$1")
what=${2:-binary}
n=${3:-1}
here=$(cd "$(dirname "$0")" && pwd)

name=$(basename "$tarball" .tar.gz)          # PBeFS-0.1.0
version=${name#PBeFS-}
out=$PWD/build-deb
rm -rf "$out"
mkdir -p "$out"
cp "$tarball" "$out/pbefs_$version.orig.tar.gz"
tar -xzf "$tarball" -C "$out"
src=$out/$name
rm -rf "$src/debian"
cp -R "$here/debian" "$src/debian"

cd "$src"
if [ "$what" = binary ]; then
  dpkg-buildpackage -b -us -uc
  echo "built: $(ls "$out"/*.deb)"
  exit 0
fi

# a PPA build: the version and the distribution of the series, the orig
# tarball only with the first upload of a version (n = 1)
series=$what
sed -i "1s/^pbefs ([^)]*) [a-z]*;/pbefs ($version-1~$series$n) $series;/" \
    debian/changelog
if [ "$n" = 1 ]; then sa=-sa; else sa=-sd; fi
dpkg-buildpackage -S $sa -d
changes=$out/pbefs_$version-1~$series${n}_source.changes
echo "to upload: dput ppa:<team>/<ppa> $changes"

###############################################################################
############################### End of deb.sh #################################
###############################################################################
