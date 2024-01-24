use cmake::Config;

fn main() {
    println!("cargo:rerun-if-env-changed=CXXFLAGS");
    println!("cargo:rerun-if-changed=../../../../icicle");

    // Base config
    let mut config = Config::new("../../../../icicle");
    config
        .define("BUILD_TESTS", "OFF")
        .define("CURVE", "bls12_377")
        .define("CMAKE_BUILD_TYPE", "Release");

    // Optional Features
    #[cfg(feature = "g2")]
    config.define("G2_DEFINED", "ON");

    // Build
    let out_dir = config
        .build_target("icicle")
        .build();

    println!("cargo:rustc-link-search={}/build", out_dir.display());
    println!("cargo:rustc-link-lib=ingo_bls12_377");

    if cfg!(feature = "bw6-761") {
        // Base config
        let mut config = Config::new("../../../../icicle");
        config
            .define("BUILD_TESTS", "OFF")
            .define("CURVE", "bw6_761")
            .define("CMAKE_BUILD_TYPE", "Release");

        // Optional Features
        #[cfg(feature = "bw6-761-g2")]
        config.define("G2_DEFINED", "ON");

        // Build
        let out_dir = config
            .build_target("icicle")
            .build();

        println!("cargo:rustc-link-search={}/build", out_dir.display());
        println!("cargo:rustc-link-lib=ingo_bw6_761");
    }

    println!("cargo:rustc-link-lib=stdc++");
    println!("cargo:rustc-link-lib=cudart");
}
