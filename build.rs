fn main() {
    println!("cargo::rerun-if-changed=src/xsf/c_api.hpp");
    println!("cargo::rerun-if-changed=src/xsf/shim.cpp");
    println!("cargo::rerun-if-changed=src/xsf/vendored");

    let out_dir = PathBuf::from(env::var_os("OUT_DIR").unwrap());

    let prefix = format!(
        "{}_{}_xsf_",
        env::var("CARGO_PKG_NAME").unwrap(),
        env::var("CARGO_PKG_VERSION").unwrap().replace(".", "_")
    );
    let prefix_flag = format!("-DPREFIX={prefix}");

    cc::Build::new()
        .cpp(true)
        .std("c++17")
        .file("src/xsf/shim.cpp")
        .flag("-Wno-unused-parameter")
        .flag(&prefix_flag)
        .compile("xsf");

    let bindings = bindgen::Builder::default()
        .header("src/xsf/c_api.hpp")
        .clang_arg(&prefix_flag)
        .allowlist_file("src/xsf/c_api.hpp")
        .blocklist_type("std::complex")
        .blocklist_type("std::complex_value_type")
        .blocklist_type("Slice")
        .generate()
        .unwrap()
        .to_string();
    let mut bindings = bindings.as_bytes();

    let mut res = Vec::new();
    while let Some((&first, rest)) = bindings.split_first() {
        if let Some(rest) = bindings.strip_prefix(b"pub fn ") {
            let (name, rest) = rest.split_at(rest.iter().position(|&b| b == b'(').unwrap());
            let name = str::from_utf8(name).unwrap();
            let stripped_name = name.strip_prefix(&prefix).unwrap();
            write!(res, "#[link_name = {name:?}] pub fn {stripped_name}").unwrap();
            bindings = rest;
            continue;
        }
        if let Some(rest) = bindings.strip_prefix(b"xsf_") {
            bindings = rest;
            continue;
        }
        res.push(first);
        bindings = rest;
    }
    fs::write(out_dir.join("xsf_bindings.rs"), &res).unwrap();
}

use std::env;
use std::fs;
use std::io::Write as _;
use std::path::PathBuf;
