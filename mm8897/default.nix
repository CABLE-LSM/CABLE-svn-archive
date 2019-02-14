let
  pkgs = import <nixpkgs> {};
  mkDerivation = import ./autotools.nix pkgs;
in mkDerivation {
  name ="cable";
  ##src =./cable2.0-trunk.tgz;
  #ctn=./CABLE-2.3.4_nix;
  ctn="CABLE-2.0_nix";
  buildInputs = with pkgs; [ mksh gfortran netcdffortran openmpi gdb];
  ncd= with pkgs; netcdffortran;
}
