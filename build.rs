use std::path::Path;

extern crate cmake;

fn main() {
    let mut config = cmake::Config::new("QuEST");

    config
        .no_build_target(true)
        .define("MULTITHREADED", "ON")
        .define("GPUACCELERATED", "OFF")
        .define("DISTRIBUTED", "OFF");

    config.define("PRECISION", "2");

    #[cfg(feature = "f32")]
    config.define("PRECISION", "1");

    #[cfg(feature = "gpu")]
    config
        .define("MULTITHREADED", "OFF")
        .define("GPUACCELERATED", "ON")
        .define("DISTRIBUTED", "OFF");

    #[cfg(feature = "mpi")]
    config
        .define("MULTITHREADED", "ON")
        .define("GPUACCELERATED", "OFF")
        .define("DISTRIBUTED", "ON");

    let dst = config.build();

    println!(
        "cargo:rustc-link-search=native={}/build/QuEST",
        dst.display()
    );
    println!("cargo:rustc-link-lib=dylib=QuEST");

    // To be able to run documentation tests, we need to work around a known
    // issue with `cargo`: [#8531](https://github.com/rust-lang/cargo/issues/8531).
    //
    // Create a link to libQuEST
    let out_dir = std::env::var_os("OUT_DIR").unwrap();
    let libfile = Path::new(&out_dir).join("build/QuEST/libQuEST.so");
    let linkfile = Path::new(&out_dir).join("../../../deps/libQuEST.so");
    let _ = std::fs::remove_file(&linkfile);
    std::os::unix::fs::symlink(libfile, linkfile).unwrap();
}
