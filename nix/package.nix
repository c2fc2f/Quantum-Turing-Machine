{
  version,
  lib,
  installShellFiles,
  rustPlatform,
  buildFeatures ? [ ],
}:

rustPlatform.buildRustPackage rec {
  pname = "qtm";

  src = lib.fileset.toSource {
    root = ../.;
    fileset = lib.fileset.unions [
      ../src
      ../build.rs
      ../Cargo.lock
      ../Cargo.toml
    ];
  };

  inherit buildFeatures;
  inherit version;

  # inject version from nix into the build
  env.NIX_RELEASE_VERSION = version;

  cargoLock.lockFile = ../Cargo.lock;

  nativeBuildInputs = [
    installShellFiles

    rustPlatform.bindgenHook
  ];

  buildInputs = [ ];

  postInstall = ''
    installShellCompletion --cmd ${meta.mainProgram} \
      --bash <($out/bin/${meta.mainProgram} completion bash) \
      --fish <($out/bin/${meta.mainProgram} completion fish) \
      --zsh <($out/bin/${meta.mainProgram} completion zsh)

    installManPage $(find target -type f -path "*/build/${pname}-*/out/*.1")
  '';

  meta = with lib; {
    description = "Binary to simulates Quantum Turing Machines (QTMs) and converts them to equivalent quantum circuits";
    mainProgram = "qtm";
    homepage = "https://github.com/c2fc2f/Quantum-Turing-Machine";
    license = licenses.mit;
    maintainers = [ maintainers.c2fc2f ];
  };
}
