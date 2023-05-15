
import click
from cloudpathlib import GSPath
from cyvcf2 import VCF, Writer
import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, output_path


vcf = VCF("~/Downloads/str_expansionhunter_v3_v3_CPG10488_eh.vcf")

fname = "out.vcf"
w = Writer(fname, vcf)

for v in vcf: 
    if v.FILTER == "PASS":
        w.write_record(v)
w.close()
vcf.close()