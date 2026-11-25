# The x86_64 Linux release bundle (tools/release_build.sh), with bonefish, in a small image:
#   docker build --build-arg VERSION=0.3.0 -t lrpipe .   (with dist/lrpipe-0.3.0-x86_64-linux/)
#   docker run --rm -v "$PWD:/data" lrpipe sv --input /data/sample.bam --vcf /data/sample.sv.vcf
#   docker run --rm -v /ref:/ref -v "$PWD:/data" lrpipe run --config /ref/lrpipe/site.toml \
#     --reads /data/sample.fastq.gz --sample S --platform ont --output /data/S
FROM debian:bookworm-slim
ARG VERSION
LABEL org.opencontainers.image.source="https://github.com/DylanCosto/lrpipe" \
      org.opencontainers.image.description="Long reads from FASTQ to VCF: alignment, SNVs and indels, SVs, copy number, tandem repeats, depth and methylation" \
      org.opencontainers.image.licenses="MIT AND BSD-3-Clause"
COPY dist/lrpipe-${VERSION}-x86_64-linux /opt/lrpipe
ENV PATH=/opt/lrpipe/bin:$PATH
ENTRYPOINT ["lrpipe"]
