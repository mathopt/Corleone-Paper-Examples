# Julia 1.12.7 with the instantiated and precompiled environment of the repository.
#   docker build -t corleone-paper .
#   docker run --rm -v "$PWD/results":/work/results corleone-paper
FROM ubuntu:24.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends curl ca-certificates \
    && rm -rf /var/lib/apt/lists/*
RUN curl -fsSL -o /tmp/julia.tar.gz \
        https://julialang-s3.julialang.org/bin/linux/x64/1.12/julia-1.12.7-linux-x86_64.tar.gz \
    && echo "4e7e9e776634d24835250de67cde39b0d4af15bc432eb20697e6be6c28ea69e8  /tmp/julia.tar.gz" | sha256sum -c - \
    && tar -xzf /tmp/julia.tar.gz -C /opt && mv /opt/julia-1.12.7 /opt/julia \
    && rm /tmp/julia.tar.gz
ENV PATH=/opt/julia/bin:$PATH JULIA_DEPOT_PATH=/opt/julia-depot
WORKDIR /work
COPY Project.toml Manifest.toml ./
RUN julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
COPY . .
CMD ["bash", "scripts/run.sh"]
