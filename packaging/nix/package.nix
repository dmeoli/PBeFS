# The derivation of PBeFS for nixpkgs, as pkgs/by-name/pb/pbefs/package.nix;
# FastFlow, which no package provides, and MPI are left out
{
  lib,
  stdenv,
  fetchurl,
  cmake,
  llvmPackages,
}:

stdenv.mkDerivation (finalAttrs: {
  pname = "pbefs";
  version = "0.1.0";

  src = fetchurl {
    url = "https://github.com/dmeoli/PBeFS/releases/download/v${finalAttrs.version}/PBeFS-${finalAttrs.version}.tar.gz";
    hash = "sha256-AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA=";
  };

  nativeBuildInputs = [ cmake ];

  buildInputs = lib.optionals stdenv.cc.isClang [ llvmPackages.openmp ];

  cmakeFlags = [
    (lib.cmakeBool "BUILD_SHARED_LIBS" (!stdenv.hostPlatform.isStatic))
    (lib.cmakeBool "PBEFS_WITH_FASTFLOW" false)
    (lib.cmakeBool "PBEFS_WITH_MPI" false)
    (lib.cmakeBool "PBEFS_BUILD_EXAMPLES" false)
    (lib.cmakeBool "PBEFS_BUILD_TOOLS" false)
    (lib.cmakeBool "PBEFS_BUILD_TESTS" finalAttrs.doCheck)
  ];

  doCheck = true;

  meta = {
    description = "Parallel best-first search in C++20: A* and branch-and-bound";
    homepage = "https://github.com/dmeoli/PBeFS";
    changelog = "https://github.com/dmeoli/PBeFS/blob/v${finalAttrs.version}/CHANGELOG.md";
    license = lib.licenses.lgpl3Only;
    maintainers = [ ];
    platforms = lib.platforms.unix;
  };
})
