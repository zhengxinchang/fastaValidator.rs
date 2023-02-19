build:
	cargo build --release --target x86_64-unknown-linux-musl

test:
	@echo "resources/test/0000231_seq.fsa"
	time ./target/x86_64-unknown-linux-musl/release/checkFastaCOVID19 resources/test/0000231_seq.fsa

	@echo "resources/test/GCF_000001405.40_GRCh38.p14_genomic.fna"
	time ./target/x86_64-unknown-linux-musl/release/checkFastaCOVID19 resources/test/GCF_000001405.40_GRCh38.p14_genomic.fna > /dev/null

	@echo "resources/test/0000231_seq.fsa.gz"
	time ./target/x86_64-unknown-linux-musl/release/checkFastaCOVID19 resources/test/0000231_seq.fsa.gz

crlf:
	@echo "resources/test/0000231_seq.fsa CRLF"
	time ./target/x86_64-unknown-linux-musl/release/checkFastaCOVID19 resources/test/0000231_seq_CRLF.fsa

h:
	time ./target/x86_64-unknown-linux-musl/release/checkFastaCOVID19 -h

ver:
	time ./target/x86_64-unknown-linux-musl/release/checkFastaCOVID19 -V


# 防止与命令相同的同名文件
.PHONY: build test 