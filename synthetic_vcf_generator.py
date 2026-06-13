import random
import argparse


def generate_af():
    r = random.random()

    if r < 0.70:
        return random.uniform(0.001, 0.05)

    elif r < 0.95:
        return random.uniform(0.05, 0.30)

    else:
        return random.uniform(0.30, 0.50)


def generate_gt(af, missing_rate=0.02):
    if random.random() < missing_rate:
        return "./."

    p00 = (1 - af) ** 2
    p01 = 2 * af * (1 - af)
    p11 = af ** 2

    r = random.random()

    if r < p00:
        return "0/0"

    elif r < p00 + p01:
        return "0/1"

    else:
        return "1/1"


def generate_pl(gt):
    if gt == "./.":
        return "."

    if gt == "0/0":
        vals = [
            0,
            random.randint(20, 80),
            random.randint(60, 200)
        ]

    elif gt == "0/1":
        vals = [
            random.randint(20, 80),
            0,
            random.randint(20, 80)
        ]

    else:
        vals = [
            random.randint(60, 200),
            random.randint(20, 80),
            0
        ]

    return ",".join(map(str, vals))


def convert_vcf(input_vcf, output_vcf, n_samples):
    with open(input_vcf) as fin, open(output_vcf, "w") as fout:

        samples = [f"SAMPLE_{i+1}" for i in range(n_samples)]

        fout.write("##fileformat=VCFv4.2\n")

        fout.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(samples)
            + "\n"
        )

        for line in fin:

            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")

            chrom = fields[0]
            pos = fields[1]
            vid = fields[2]
            ref = fields[3]
            alt = fields[4]

            qual = fields[5] if len(fields) > 5 else "."
            filt = fields[6] if len(fields) > 6 else "PASS"

            af = generate_af()

            genotypes = []

            for _ in range(n_samples):

                gt = generate_gt(af)
                pl = generate_pl(gt)

                if gt == "./.":
                    genotypes.append("./.:.")
                else:
                    genotypes.append(f"{gt}:{pl}")

            row = [
                chrom,
                pos,
                vid,
                ref,
                alt,
                qual,
                filt,
                f"AF={af:.4f}",
                "GT:PL"
            ] + genotypes

            fout.write("\t".join(row) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="GWAS VCF-like input"
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Synthetic VCF output"
    )

    parser.add_argument(
        "-n",
        "--samples",
        type=int,
        default=100,
        help="Number of synthetic samples"
    )

    args = parser.parse_args()

    convert_vcf(
        args.input,
        args.output,
        args.samples
    )
