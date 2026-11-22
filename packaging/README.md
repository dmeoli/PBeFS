# Packaging

This directory holds the recipes with which PBeFS is distributed by the package
managers, one subdirectory per channel, together with the scripts that fill
them at each release and a small project that uses the installed library
(`test/`), which the recipes build against the package they have just made. All
the recipes take the same source, i.e., the tarball `PBeFS-X.Y.Z.tar.gz` that
the workflow `release.yml` attaches, together with its sha256, to the GitHub
release of the tag `vX.Y.Z` (cf. the section on the release below). Git fills
in the `VERSION.txt` of the tarball with the tag, and CMake takes the version
of the package from it (see `cmake/DeriveVersion.cmake`), also when the tarball
is unpacked inside another git repository, as in the build trees of vcpkg.

Since PBeFS is a compiled library, its packages depend on the architecture.
Each of them installs the headers under `include/PBeFS`, the library `PBeFS`
(static or shared, as the channel prefers), the CMake package `PBeFS` (target
`PBeFS::PBeFS`) under `lib/cmake/PBeFS` and the pkg-config file `pbefs.pc`
under `lib/pkgconfig`. Each also depends on the package of ADS of its channel
(`ads`, `cpp-ads` on conda-forge, `libads-dev` on Debian), whose headers those
of PBeFS include. In the following table, "Channel" is the package manager,
"Files" the files of its recipe in this directory, and "Name" the name of the
package on that channel.

| Channel        | Files                              | Name                |
|----------------|------------------------------------|---------------------|
| vcpkg          | `vcpkg/ports/pbefs/`               | `pbefs`             |
| conda-forge    | `conda/recipe/recipe.yaml`         | `pbefs`             |
| Homebrew       | `homebrew/Formula/pbefs.rb`        | `pbefs`             |
| Debian, PPA    | `debian/`, `deb.sh`                | `libpbefs-dev`      |
| Conan Center   | `conan/`                           | `pbefs`             |
| Spack          | `spack/packages/pbefs/`            | `pbefs`             |
| Nix            | `nix/package.nix`, `flake.nix`     | `pbefs`             |
| AUR            | `aur/PKGBUILD`, `aur/.SRCINFO`     | `pbefs`             |
| MSYS2          | `msys2/PKGBUILD`                   | `mingw-w64-*-pbefs` |
| xmake-repo     | `xmake/packages/p/pbefs/`          | `pbefs`             |
| Meson WrapDB   | `wrapdb/`                          | `pbefs`             |

As far as their indexes showed when the recipes were written, the name `pbefs`
was free on all these channels (and `libpbefs-dev` in Debian and Ubuntu). In
each recipe the parallel runtimes stay optional, and none of them is found by
chance on the machine that builds the package. OpenMP, which comes with GCC and
MSVC, gives the target `PBeFS::OpenMP` where the compiler has it (vcpkg and
Spack take it as a feature and as a variant), while MPI gives `PBeFS::MPI` only
with the feature `mpi` of vcpkg and with the variant `+mpi` of Spack. FastFlow
and Sylvan (and with the latter the library of the symbolic searches,
`PBeFS::Symbolic`) are left out of all of them, since no channel packages them
and the build would download them at configure time, which the build farms of
the channels forbid. The headers of all the searches are installed anyway, and
a project that has MPI can use `HDAStarMPI` by linking `PBeFS::PBeFS` together
with `MPI::MPI_CXX`; in fact, the compiled part of the library depends on no
runtime. How the installed CMake package imports the target of a runtime is
described in the section on the use of the library of the README of the
repository.

On each push that touches the packaging, the build or the sources of the
library, the workflow `packaging.yml` builds and tests the Homebrew formula on
macOS, the vcpkg port and the MSYS2 package on Windows, the AUR package on Arch
Linux, the Debian package on Ubuntu and the port of WrapDB, each from the
tarball of the commit.

The Bazel Central Registry is left out. A module of PBeFS there would give only
the library on the C++ threads, since the rules of Bazel have no portable way
to OpenMP and MPI, and it would need a `MODULE.bazel` with a version written by
hand (or an overlay in the registry), whereas the version of PBeFS comes from
its tag. WrapDB, on the contrary, keeps its port as an overlay in its
repository, and that port is a single `meson.build` (cf. the section on
WrapDB).

## The release

In the recipes, the tarball of the release has a placeholder checksum, which is
filled once the tarball exists, as follows.

1. The section `[Unreleased]` of `CHANGELOG.md` becomes `## [X.Y.Z] - <date>`,
   and the commit is tagged `vX.Y.Z` and pushed. Then `release.yml` checks that
   the `VERSION.txt` of the tarball holds the tag, and publishes the release,
   with the notes taken from that section, with `PBeFS-X.Y.Z.tar.gz` and
   `PBeFS-X.Y.Z.tar.gz.sha256` attached.
2. `packaging/bump.sh vX.Y.Z` downloads the tarball, checks its `VERSION.txt`
   again, and writes its version and its checksums into the recipes (the
   sha256, the sha512 of vcpkg and the SRI hash of Nix). Those that keep all
   the versions (Conan, Spack, xmake) get a new entry, and the others move to
   the new version with their release number back to 1; `debian/changelog` gets
   a new entry, the wrap of WrapDB a new release `X.Y.Z-1`, and the snippets of
   FetchContent and CPM.cmake in the README of the repository take the new
   version. This result is committed on `main`.
3. Each channel below is then updated from the committed recipes. Only the
   first release needs a submission to each channel, since several of them
   (conda-forge, Homebrew, Spack, Nix, vcpkg) have bots or reviewers that
   follow the later releases, possibly with a pull request to approve.

Note that the tarball must not change after the release, since all the
checksums depend on it; hence, a release that has to be redone takes a new tag.

## vcpkg

The port is tested as an overlay with

    vcpkg install "pbefs[openmp]" --overlay-ports=packaging/vcpkg/ports

and the submission to the curated registry goes through a fork of
`microsoft/vcpkg`,

    git clone https://github.com/<user>/vcpkg && cd vcpkg
    ./bootstrap-vcpkg.sh
    cp -R <PBeFS>/packaging/vcpkg/ports/pbefs ports/
    ./vcpkg format-manifest ports/pbefs/vcpkg.json
    ./vcpkg install "pbefs[openmp,mpi]"
    git add ports/pbefs && git commit -m "[pbefs] new port"
    ./vcpkg x-add-version pbefs
    git add versions && git commit --amend --no-edit

followed by a pull request to `microsoft/vcpkg` with the checklist of new ports
filled in. Its CI builds the port on all its triplets, static and dynamic (the
DLL exports all the symbols of the library, by `WINDOWS_EXPORT_ALL_SYMBOLS`),
and the review usually takes one to three weeks. Reviewers may request
`vcpkg_from_github` (the archive of the tag, with `REF v${VERSION}` and the
SHA512 of that archive) in place of the asset of the release, which is a
two-line change of `portfile.cmake`. GitHub makes that archive with git archive
too, but whether its `VERSION.txt` holds the tag, which CMake needs, is to be
checked then. A later release is a pull request with the new `portfile.cmake`,
`vcpkg.json` and `x-add-version`.

## conda-forge

A new package enters through `conda-forge/staged-recipes`,

    git clone https://github.com/<user>/staged-recipes
    mkdir staged-recipes/recipes/pbefs
    cp packaging/conda/recipe/recipe.yaml staged-recipes/recipes/pbefs/

with a pull request, whose CI builds the recipe with rattler-build on Linux,
macOS and Windows and runs its tests (the content of the package, then
`packaging/test` with CMake and pkg-config); the linter is
`conda-smithy recipe-lint --conda-forge`. This recipe builds the shared
library, which conda-forge prefers to the static one. One then requests a
review with a comment `@conda-forge/help-c-cpp, ready for review!`, which
typically takes one to four weeks. After the merge, the feedstock
`conda-forge/pbefs-feedstock` is created with `dmeoli` as its maintainer, and
the autotick bot, which watches the releases on GitHub, opens the pull request
of each later release.

## Homebrew

`homebrew-core` accepts only projects with some notability (at least 30 forks,
30 watchers or 75 stars on GitHub), and it usually refuses the new ones
submitted by their author. Hence, the formula goes first to a separate tap, a
GitHub repository `dmeoli/homebrew-pbefs` with the formula in
`Formula/pbefs.rb`, after which

    brew install dmeoli/pbefs/pbefs

installs the shared library on macOS and on Linux (with `PBeFS::OpenMP` on
Linux, where the compiler of Homebrew is GCC, and without it on macOS, whose
Apple clang has no OpenMP). Before pushing the formula one runs
`brew style Formula/pbefs.rb`,
`brew install --build-from-source dmeoli/pbefs/pbefs`,
`brew test dmeoli/pbefs/pbefs` and
`brew audit --strict --online dmeoli/pbefs/pbefs`; a workflow running
`brew test-bot` in the tap can also build the bottles. Once PBeFS is notable,
the same formula goes to `homebrew-core` as `Formula/p/pbefs.rb` with a pull
request, whose CI builds it on the supported versions of macOS, and later
versions are then bumped by `brew bump-formula-pr`.

## Debian and Ubuntu

Here the Debian sources stay in `packaging/debian`, outside the root of the
tarball, as Debian requires of upstream projects, and `deb.sh` puts them into
the unpacked tarball. Its package `libpbefs-dev` holds the headers, the static
library and the CMake and pkg-config files; this spares a package of the shared
library with its symbols while the binary interface changes at each minor
release. A local build, with the tests of the library run by `dh_auto_test`
(skipped with `DEB_BUILD_OPTIONS=nocheck`), is

    packaging/deb.sh PBeFS-X.Y.Z.tar.gz binary
    lintian -I --pedantic build-deb/*.changes

For a Launchpad PPA (say `ppa:dmeoli/pbefs`), once it exists under the account
with the GPG key registered, the signed source package of a series and its
upload are

    packaging/deb.sh PBeFS-X.Y.Z.tar.gz resolute
    dput ppa:dmeoli/pbefs build-deb/pbefs_X.Y.Z-1~resolute1_source.changes

and the same with `noble` for Ubuntu 24.04, whose GCC 13 compiles the library
and the tests. A second upload of the same version to the same series takes
`n = 2` (`deb.sh <tarball> <series> 2`), which leaves the original tarball out.
Launchpad builds the package within an hour and mails the result. Entering
Debian itself takes an Intent To Package (a bug against `wnpp`, filed with
`reportbug wnpp`), whose number goes into the first entry of `debian/changelog`
as `Closes: #<number>`, an upload to `mentors.debian.net` and a Debian
Developer who sponsors it; the review there is thorough and may take months.

## Conan Center

The recipe follows the layout of `conan-io/conan-center-index`, where it goes
as `recipes/pbefs/` (the files `config.yml` and `all/` of `packaging/conan`),
and

    conan create recipes/pbefs/all --version X.Y.Z --build=missing

builds the package and runs `test_package` with a profile in C++20 (with C++17
the recipe declares the configuration invalid). Since the files of CMake that
Conan gives to a consumer are those of its generator `CMakeDeps` instead of the
package of the project, the recipe gives the library (`PBeFS::PBeFS`) built
without the runtimes. Its pull request is built by the CI of Conan Center on
Linux, macOS and Windows, and the review of a new recipe is slow (several weeks
is common). A later version adds an entry to `config.yml` and
`all/conandata.yml`, which `bump.sh` writes.

## Spack

The package goes to `spack/spack-packages` as
`repos/spack_repo/builtin/packages/pbefs/package.py`. In a clone of Spack with
that repository,

    spack style --fix && spack audit packages pbefs
    spack install --test=root pbefs +openmp +mpi

checks the recipe, then builds the library and runs its tests (with `+mpi`,
also the searches over MPI). Its pull request is reviewed by the maintainers of
Spack, usually within a few weeks, and `maintainers("dmeoli")` makes the user
the reviewer of the later changes to the package. Each new release adds a line
`version(...)`, which `bump.sh` writes.

## Nix

The flake at the root of the repository builds the checkout, and hence
`nix build github:dmeoli/PBeFS` and `nix flake check github:dmeoli/PBeFS` work
with no submission at all; since the checkout has no tag, the version of its
package is 0.0.0. For nixpkgs, in a fork of `NixOS/nixpkgs`,

    mkdir -p pkgs/by-name/pb/pbefs
    cp <PBeFS>/packaging/nix/package.nix pkgs/by-name/pb/pbefs/
    nixfmt pkgs/by-name/pb/pbefs/package.nix
    nix-build -A pbefs

with the user added to `maintainers/maintainer-list.nix` in a first commit
("maintainers: add dmeoli") and listed in `meta.maintainers`, and the package
in a second commit titled "pbefs: init at X.Y.Z". In nixpkgs, the review of a
new package may take weeks, and later releases are then bumped by `r-ryantm`.

## AUR

The AUR has no review. With an account on `aur.archlinux.org` and an SSH key
registered there,

    git clone ssh://aur@aur.archlinux.org/pbefs.git aur-pbefs
    cp packaging/aur/PKGBUILD aur-pbefs/ && cd aur-pbefs
    makepkg --printsrcinfo > .SRCINFO
    makepkg -si && namcap PKGBUILD pbefs-*.pkg.tar.zst
    git add PKGBUILD .SRCINFO && git commit -m "pbefs X.Y.Z" && git push

publishes the package at once (the `.SRCINFO` of `bump.sh` is the same file,
but on Arch `makepkg --printsrcinfo` is the reference). This package installs
the shared library, as Arch prefers, and its `check()` builds and runs the
tests of the library, which `makepkg --nocheck` skips.

## MSYS2

The package goes to `msys2/MINGW-packages` as `mingw-w64-pbefs/PKGBUILD`. In an
MSYS2 shell,

    cd mingw-w64-pbefs && makepkg-mingw --cleanbuild --syncdeps --force

builds it for the current environment (UCRT64, CLANG64, and so on). In turn,
the CI of MINGW-packages builds the environments listed in `mingw_arch`, and a
new package is usually merged within days.

## xmake-repo

The package goes to `xmake-io/xmake-repo`, on the branch `dev`, as
`packages/p/pbefs/xmake.lua`. In a clone of that repository,

    xmake l scripts/test.lua -vD --shallow pbefs

builds it and runs its `on_test`, which compiles a snippet against
`PBeFS/PBeFS.h` in C++20. Its pull request goes against `dev`, and the
maintainers usually merge it within days; a later version is a new line
`add_versions(...)`, which `bump.sh` writes.

## Meson WrapDB

WrapDB takes a project that does not build with Meson through a port kept in
its repository, `mesonbuild/wrapdb`, as an overlay on the tarball, and
`packaging/wrapdb` holds the files of the port with the layout of that
repository. These are the wrap `subprojects/pbefs.wrap` (the tarball, its
sha256, the overlay and the dependency `pbefs` it provides), the overlay
`subprojects/packagefiles/pbefs/meson.build`, and the entry of `releases.json`,
whose versions are `X.Y.Z-N`, N being the revision of the port. This overlay
compiles the library from the four sources of `src` with the threads (static
with MSVC, which exports no symbol unless told to), and gives the dependency
`pbefs` with the headers and the threads, the installation of the headers and
of `pbefs.pc`, and `packaging/test/main.cpp` as its test. Note that the
runtimes of OpenMP, MPI and FastFlow are headers, which a project takes
together with its dependency on the runtime. In a clone of a fork of WrapDB,

    cp <PBeFS>/packaging/wrapdb/subprojects/pbefs.wrap subprojects/
    cp -R <PBeFS>/packaging/wrapdb/subprojects/packagefiles/pbefs \
          subprojects/packagefiles/
    # the entry "pbefs" of <PBeFS>/packaging/wrapdb/releases.json goes
    # into releases.json, in the order of the names
    ./tools/format.py
    git add -A && git commit -m "pbefs: add X.Y.Z"
    ./tools/sanity_checks.py
    meson setup builddir -Dwraps=pbefs && meson test -C builddir

with a pull request to `mesonbuild/wrapdb` after them, whose CI repeats the
checks and the build on Linux, macOS and Windows before a maintainer merges it.
After the merge, a project takes PBeFS as shown in the README of the
repository. A later version is a pull request with the new wrap and the new
entry of `releases.json`, which `bump.sh` writes, and a fix of the port only
takes a new revision, `X.Y.Z-2`.

## CMake without a package manager

A CMake project can also take PBeFS from its sources at configure time, with
FetchContent or CPM.cmake, as shown in the section on its use in the README of
the repository.
