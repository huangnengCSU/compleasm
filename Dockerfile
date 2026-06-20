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

# Install Python dependencies and sepp
RUN pip3 install --break-system-packages pandas sepp

# Install compleasm from GitHub
RUN pip3 install --break-system-packages \
    git+https://github.com/huangnengCSU/compleasm.git@0.2.6

ENTRYPOINT ["compleasm"]
