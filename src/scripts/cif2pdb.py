# -*- coding: utf-8 -*-
"""mmCIF -> PDB 转换 (BioPython), 用于对接受体准备"""
import sys
from Bio.PDB import MMCIFParser, PDBIO

cif, pdb = sys.argv[1], sys.argv[2]
s = MMCIFParser(QUIET=True).get_structure("x", cif)
io = PDBIO()
io.set_structure(s)
io.save(pdb)
import os
print("saved", pdb, os.path.getsize(pdb), "bytes")
