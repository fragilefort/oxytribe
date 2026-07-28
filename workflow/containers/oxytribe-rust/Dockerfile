FROM rust:1-slim-bookworm AS builder

WORKDIR /usr/src/app

COPY Cargo.toml Cargo.lock ./
COPY crates ./crates

RUN cargo build --release

FROM debian:bookworm-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/src/app/target/release/combine_edit_sites /usr/local/bin/
COPY --from=builder /usr/src/app/target/release/find_edit_sites /usr/local/bin/
COPY --from=builder /usr/src/app/target/release/sam_to_map /usr/local/bin/

ENV PATH="/usr/local/bin:${PATH}"
