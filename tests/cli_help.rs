use assert_cmd::Command;

#[test]
fn cli_help_smoke() {
    let mut cmd = Command::cargo_bin("kira-proteoqc").unwrap();
    cmd.arg("--help");
    cmd.assert().success();
}
