extern crate cmake;

fn main() {
    let dst = cmake::Config::new("QuEST")
        .no_build_target(true)
        .define("MULTITHREADED", "ON")
        .define("GPUACCELERATED", "OFF")
        .define("DISTRIBUTED", "OFF")
        .build();

    println!(
        "cargo:rustc-link-search=native={}/build/QuEST",
        dst.display()
    );
    println!("cargo:rustc-link-lib=dylib=QuEST");
}
