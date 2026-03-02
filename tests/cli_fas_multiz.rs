#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::fs;
use tempfile::TempDir;

#[test]
fn command_fas_multiz_core() {
    let tempdir = TempDir::new().unwrap();
    let out_path = tempdir.path().join("merged.fas");
    let out_str = out_path.to_str().unwrap();

    PgrCmd::new()
        .args(&[
            "fas",
            "multiz",
            "-r",
            "S288c",
            "tests/fas/S288cvsRM11_1a.slice.fas",
            "tests/fas/S288cvsYJM789.slice.fas",
            "tests/fas/S288cvsSpar.slice.fas",
            "-o",
            out_str,
        ])
        .assert()
        .success();

    assert!(out_path.is_file());
    let content = fs::read_to_string(out_path).unwrap();
    assert!(content.lines().count() > 0);

    tempdir.close().unwrap();
}

#[test]
fn command_fas_multiz_affine_gap() {
    let tempdir = TempDir::new().unwrap();
    let out_path = tempdir.path().join("merged_affine.fas");
    let out_str = out_path.to_str().unwrap();

    PgrCmd::new()
        .args(&[
            "fas",
            "multiz",
            "-r",
            "S288c",
            "tests/fas/S288cvsRM11_1a.slice.fas",
            "tests/fas/S288cvsYJM789.slice.fas",
            "tests/fas/S288cvsSpar.slice.fas",
            "--gap-open",
            "400",
            "--gap-extend",
            "30",
            "-o",
            out_str,
        ])
        .assert()
        .success();

    assert!(out_path.is_file());
    let content = fs::read_to_string(out_path).unwrap();
    assert!(content.lines().count() > 0);

    tempdir.close().unwrap();
}

#[test]
fn command_fas_multiz_custom_matrix() {
    let tempdir = TempDir::new().unwrap();
    let out_path = tempdir.path().join("merged_matrix.fas");
    let out_str = out_path.to_str().unwrap();

    PgrCmd::new()
        .args(&[
            "fas",
            "multiz",
            "-r",
            "S288c",
            "tests/fas/S288cvsRM11_1a.slice.fas",
            "tests/fas/S288cvsYJM789.slice.fas",
            "tests/fas/S288cvsSpar.slice.fas",
            "--score-matrix",
            "hoxd55",
            "-o",
            out_str,
        ])
        .assert()
        .success();

    assert!(out_path.is_file());
    let content = fs::read_to_string(out_path).unwrap();
    assert!(content.lines().count() > 0);

    tempdir.close().unwrap();
}
