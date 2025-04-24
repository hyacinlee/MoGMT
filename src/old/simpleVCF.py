import re
import sys

def simplify_and_filter_vcf(input_file, output_file):
    with open(input_file) as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue

            parts = line.strip().split('\t')
            format_fields = parts[8].split(':')

            try:
                gt_index = format_fields.index("GT")
                ad_index = format_fields.index("AD")
                dp_index = format_fields.index("DP")
            except ValueError:
                continue

            new_sample_fields = []
            all_genotypes = set()
            skip_line = False

            for i in range(9, len(parts)):
                fields = parts[i].split(':')
                if len(fields) <= max(gt_index, ad_index, dp_index):
                    new_sample_fields.append('./.:.:.')
                    all_genotypes.add('./.')
                    continue

                raw_gt = fields[gt_index]
                if '.' in raw_gt:
                    new_sample_fields.append('./.:.:.')
                    all_genotypes.add('./.')
                    continue

                alleles = re.split(r'[|/]', raw_gt)
                if len(alleles) != 2:
                    new_sample_fields.append('./.:.:.')
                    all_genotypes.add('./.')
                    continue

                alleles = sorted(alleles)
                new_gt = '/'.join(alleles)

                ad = fields[ad_index]
                dp = fields[dp_index]

                new_sample_fields.append(f"{new_gt}:{ad}:{dp}")
                all_genotypes.add(new_gt)

            if len(all_genotypes) <= 1:
                continue

            parts[8] = "GT:AD:DP"
            parts[9:] = new_sample_fields
            fout.write('\t'.join(parts) + '\n')


simplify_and_filter_vcf(sys.argv[1], sys.argv[2])

