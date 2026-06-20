FROM debian:bookworm-slim

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 python3-pip python3-venv \
        hmmer \
        git \
        wget \
        build-essential \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install miniprot from source
RUN wget -q https://github.com/lh3/miniprot/releases/download/v0.13/miniprot-0.13_x64-linux.tar.bz2 \
    && tar -xjf miniprot-0.13_x64-linux.tar.bz2 \
    && cp miniprot-0.13_x64-linux/miniprot /usr/local/bin/miniprot \
    && rm -rf miniprot-0.13_x64-linux miniprot-0.13_x64-linux.tar.bz2

# Install Python dependencies
RUN pip3 install --break-system-packages pandas

# Install sepp from source
RUN pip3 install --break-system-packages git+https://github.com/smirarab/sepp.git \
    && mkdir -p /opt/sepp-home \
    && DIST_INFO=$(find /usr/local/lib -name "sepp-*.dist-info" -type d) \
    && echo "/opt/sepp-home" > "$DIST_INFO/home.path"

# Install compleasm from GitHub
RUN pip3 install --break-system-packages \
    git+https://github.com/huangnengCSU/compleasm.git@0.2.6

CMD ["compleasm"]
