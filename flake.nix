# PBeFS as a Nix flake: `nix build github:dmeoli/PBeFS` installs the library,
# the headers, the CMake package and pbefs.pc; `nix flake check` also runs
# the tests. The checkout has no tag, and hence its version is 0.0.0.
{
  description = "Parallel best-first search in C++20: A* and branch-and-bound";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs =
    { self, nixpkgs }:
    let
      systems = [
        "x86_64-linux"
        "aarch64-linux"
        "x86_64-darwin"
        "aarch64-darwin"
      ];
      each =
        f: nixpkgs.lib.genAttrs systems (s: f nixpkgs.legacyPackages.${s});
    in
    {
      packages = each (pkgs: {
        default =
          (pkgs.callPackage ./packaging/nix/package.nix { }).overrideAttrs
            {
              version = "main";
              src = self;
            };
      });
      checks = each (pkgs: {
        default = self.packages.${pkgs.stdenv.hostPlatform.system}.default;
      });
    };
}
