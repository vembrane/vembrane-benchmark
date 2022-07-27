md5sums = {}
for md5sum_file in snakemake.input.md5sums:
    with open(md5sum_file, "rt") as f:
        s = f.read()
        md5sums[md5sum_file] = s.strip()
identical = len(set(md5sums.values())) == 1
with open(snakemake.output.f, "wt") as f:
    f.write(f"#{identical}\n")
    for md5sum_file, md5sum in md5sums.items():
        f.write(f"{md5sum_file}\t{md5sum}\n")