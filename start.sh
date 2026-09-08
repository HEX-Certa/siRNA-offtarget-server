cargo build --release
DATA_DIR=data LISTEN_ADDR=0.0.0.0:8080 ./target/release/sirna-offtarget-checker
