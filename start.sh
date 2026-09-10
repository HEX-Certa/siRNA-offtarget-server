cargo build --release
DATA_DIR=data LISTEN_ADDR=0.0.0.0:9100 ./target/release/sirna-offtarget-checker
