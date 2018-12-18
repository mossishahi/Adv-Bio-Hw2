from random import randint, choice, seed

from numpy.random import uniform, seed as nseed


class Sequence:
    def __init__(self, seq_str, name="", quality_str=""):
        self.seq_str = seq_str
        self.size = len(seq_str)
        self.name = name
        self.quality_str = quality_str

    @staticmethod
    def read_big_sequence(file_name):
        file = open(file_name)
        stream = file.read()
        file.close()
        reads = stream.split(">")
        result = [Sequence("")] * 0
        for read_i in range(1 if stream[0] == ">" else 0, len(reads)):
            read = reads[read_i]
            first_new_line = read.find("\n") if stream[0] == ">" else 0
            name = read[:first_new_line]
            seq_str = read[first_new_line + 1:].replace("\n", "").upper()
            result.append(Sequence(seq_str, name))
        return result


# not_Ns = list(range(10000, 44821)) + list(range(94821, 133871)) + list(range(222346, 226276)) + list(
#     range(226351, 1949345)) + list(range(2132994, 2137388)) + list(range(2137488, 37099262)) + list(
#     range(37285837, 49348394)) + list(range(49528394, 50228964)) + list(range(50278964, 58555579)) + list(
#     range(58605579, 58736186)) + list(range(58736286, 58804505)) + list(range(58804605, 58988881)) + list(
#     range(58988981, 59064833)) + list(range(59064933, 59477728)) + list(range(59477828, 59487434)) + list(
#     range(59487534, 59532975)) + list(range(59533075, 60550108)) + list(range(60550208, 61173171)) + list(
#     range(61173271, 61212316)) + list(range(61212416, 62068000)) + list(range(62068100, 62107293)) + list(
#     range(62107393, 62412542)) + list(range(62462542, 114281198)) + list(range(114331198, 115738949)) + list(
#     range(115838949, 116557779)) + list(range(116595566, 120879381)) + list(range(120929381, 144425606)) + list(
#     range(144475606, 156030895))
seed(0)
nseed(0)

chrX = Sequence.read_big_sequence("chrX.fa")[0]
change_poses = uniform(0, chrX.size - 5, 100000)

cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
vcf_file = open('Sample.vcf', 'w')
vcf_file.write("##fileformat=VCFv4.0\n")
vcf_file.write("#" + "\t".join(cols) + "\n")

change_poses_set = set()
changes = []
bases = ["A", "C", "G", "T"]
# for i in range(chrX.size - 1):
#     if chrX.seq_str[i] != "N" and chrX.seq_str[i + 1] == "N":
#         print("start of N at %d" % (i + 1))
#     if chrX.seq_str[i] == "N" and chrX.seq_str[i + 1] != "N":
#         print("end of N at %d" % i)
# print(chrX.size)
last_unused_change_pos = 0
for i in range(3):  # INSERT, SNP, DELETE
    for j in range(1000):
        change_pos = int(change_poses[last_unused_change_pos])
        while change_pos in change_poses_set or chrX.seq_str[change_pos] == 'N' or chrX.seq_str[change_pos + 5] == 'N':
            last_unused_change_pos += 1
            change_pos = int(change_poses[last_unused_change_pos])
        last_unused_change_pos += 1
        change_poses_set.add(change_pos)
        ref_bases = chrX.seq_str[change_pos:change_pos + (1 if i < 2 else randint(1, 5))]
        new_bases = "" if i == 1 else ref_bases[0] if i == 0 else "."
        if i == 1:
            bases.remove(ref_bases[0])
            new_bases = choice(bases)
            bases.append(ref_bases[0])
        if i == 0:
            new_bases += "".join([choice(bases) for _ in range(randint(1, 5))])
        vcf_file.write("\t".join(["chrX", str(change_pos + 1), ".", ref_bases, new_bases, "41", "PASS", "."]) + "\n")
        changes.append((change_pos, ref_bases, new_bases))
vcf_file.close()

changes.sort()
new_seq_str = list(chrX.seq_str)
for pos, ref, new_base in changes:
    if "".join(new_seq_str[pos:pos + len(ref)]) != ref:
        print("salam")
    for i in range(len(ref)):
        new_seq_str[pos + i] = ""
    if new_base != ".":
        new_seq_str[pos] = new_base
new_seq_str = "".join(new_seq_str)
print(len(new_seq_str), chrX.size)
f = open("Sample.fa", "w")
f.write(">" + chrX.name + "\n")
f.write(new_seq_str)
