from typing import List

class Clash:
    def __init__(self, line):
        self.chain1, self.resnum1, self.res1, self.atom1, \
            self.chain2, self.resnum2, self.res2, self.atom2, \
            self.degree = \
                line[3], line[4:8], line[10:13], line[14:18], \
                line[24], line[25:29], line[31:34], line[35:39], float(line[48:53])
    
    def __str__(self):
        return "Clash between residue {}:{} ({} - {}) and {}:{} ({} - {}) of severity {}" \
            .format(self.chain1, self.resnum1.strip(), self.res1.strip(), self.atom1.strip(), \
            self.chain2, self.resnum2.strip(), self.res2.strip(), self.atom2.strip(), \
            self.degree)

class BondLength:
    pass
class BondAngle:
    pass
class Suite:
    pass

def bond_lengths_from_file_contents(file_contents: List[str]) -> List[BondLength]:
    look = False
    for line in file_contents:
        if "----------Bond lengths----------" in line:
            look = True
        if "----------Bond angles----------" in line: break

    return []

def bond_angles_from_file_contents(file_contents: List[str]) -> List[BondAngle]: pass

def clashes_from_file_contents(file_contents: List[str]) -> List[Clash]:
    look = False
    clashes = []
    for line in file_contents:
        if "----------Bad clashes----------" in line:
            look = True
        else:
            if look and len(line) > 50: clashes.append(Clash(line))
        if "clashscore" in line: break
    return clashes

def clashscore_from_file_contents(file_contents: List[str]) -> float:
    for line in file_contents:
        if "clashscore" in line: return float(line.split()[2])
        if "Clashscore" in line: return float(line.split()[2])
    return 0

def suites_from_file_contents(file_contents: List[str]) -> List[Suite]: pass
def rwork_from_file_contents(file_contents: List[str]) -> float: pass
def rfree_from_file_contents(file_contents: List[str]) -> float: pass



class Result:
    def __init__(self, file_contents: List[str]):
        #self.n_atoms = natoms_from_file_contents(file_contents)
        #self.
        self.bond_lengths        = bond_lengths_from_file_contents(file_contents)
        self.bond_angles         = bond_angles_from_file_contents(file_contents)
        #self.dihedral_angles     = bond_angles_from_file_contents(file_contents)
        #self.chiral_volumes      = chiral_volumes_from_file_contents(file_contents)
        #self.planar_groups       = planar_groups_from_file_contents(file_contents)
        self.clashes             = clashes_from_file_contents(file_contents)
        self.clashscore          = clashscore_from_file_contents(file_contents)
        self.suites              = suites_from_file_contents(file_contents)
        #self.puckers             = puckers_from_file_contents(file_contents)
        self.rwork               = rwork_from_file_contents(file_contents)
        self.rfree               = rfree_from_file_contents(file_contents)

if __name__ == '__main__':
    result = Result(open('/Users/amw579/projects/ERRASER_Projects/4mgm/molprobity.out').readlines())
    print([str(c) for c in result.clashes])
