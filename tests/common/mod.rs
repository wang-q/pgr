use assert_cmd::assert::Assert;
use assert_cmd::Command;

pub struct PgrCmd {
    cmd: Command,
    stdin: Option<String>,
}

impl PgrCmd {
    pub fn new() -> Self {
        let mut cmd = Command::cargo_bin("pgr").unwrap();
        cmd.env("RUST_BACKTRACE", "1");
        Self { cmd, stdin: None }
    }

    #[allow(dead_code)]
    pub fn stdin<S: Into<String>>(mut self, input: S) -> Self {
        self.stdin = Some(input.into());
        self
    }

    #[allow(dead_code)]
    pub fn args(mut self, args: &[&str]) -> Self {
        self.cmd.args(args);
        self
    }

    pub fn run(mut self) -> (String, String) {
        if let Some(input) = self.stdin {
            self.cmd.write_stdin(input);
        }

        let output = self.cmd.output().expect("Failed to execute command");

        let stdout = String::from_utf8(output.stdout).expect("Invalid UTF-8 in stdout");
        let stderr = String::from_utf8(output.stderr).expect("Invalid UTF-8 in stderr");

        (stdout, stderr)
    }

    #[allow(dead_code)]
    pub fn assert(mut self) -> Assert {
        if let Some(input) = self.stdin {
            self.cmd.write_stdin(input);
        }
        self.cmd.assert()
    }

    #[allow(dead_code)]
    pub fn run_fail(self) -> (String, String) {
        let assert = self.assert().failure();
        let output = assert.get_output();
        let stdout = String::from_utf8(output.stdout.clone()).expect("Invalid UTF-8 in stdout");
        let stderr = String::from_utf8(output.stderr.clone()).expect("Invalid UTF-8 in stderr");
        (stdout, stderr)
    }
}

#[allow(dead_code)]
pub fn assert_close(actual: f64, expected: f64, epsilon: f64) {
    if (actual - expected).abs() > epsilon {
        panic!(
            "Assertion failed: {} is not close to {} (epsilon: {}, diff: {})",
            actual,
            expected,
            epsilon,
            (actual - expected).abs()
        );
    }
}
