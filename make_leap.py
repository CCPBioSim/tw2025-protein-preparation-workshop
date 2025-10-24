#! /usr/bin/env python3
'''
A utility to help generate tleap scripts for AMBER parameterization.

Copyright 2025 Charlie Laughton'''

from argparse import ArgumentParser
import sys
import mdtraj as mdt
from pathlib import Path
from crossflow.tasks import SubprocessTask
from crossflow.filehandling import FileHandler
import shutil

__version__ = '1.0.0'


def _aliased(cmd):
    '''
    Little utility to see if a comand is aliased
    '''
    alias = SubprocessTask('alias {cmd}')
    alias.set_inputs(['cmd'])
    alias.set_outputs(['STDOUT'])
    al = alias(cmd)
    return al


def _check_available(cmd):
    '''
    Little utility to check a required command is available

    '''
    if shutil.which(cmd) is None:
        if _aliased(cmd) == '':
            raise FileNotFoundError(f'Error: cannot find the {cmd} command')


def _check_exists(filename):
    '''
    Little utility to check if a required file is present

    '''
    if not Path(filename).exists():
        raise FileNotFoundError(f'Error: cannot find required file {filename}')


def leap(amberpdb, ff, het_names=None, solvate=None, padding=10.0, het_dir='.',
         n_na=0, n_cl=0, script_only=False):
    '''
    Parameterize a molecular system using tleap (or just generate the script).

    Args:
       amberpdb (str): An Amber-compliant PDB file
       ff (list): The force fields to use.
       het_names (list): List of parameterised heterogens
       solvate (str or None): type of periodic box ('box', 'cube', or 'oct')
       padding (float): Clearance between solute and any box edge (Angstroms)
       het_dir (str or Path): location of the directory containing heterogen
                              parameters
       n_na (int): number of Na+ ions to add (0 = minimal salt)
       n_cl (int): number of Cl- ions to add (0 = minimal salt)
       script_only (bool): if True, only generate the tleap script
                            and do not run tleap

    Returns:
         (FileHandle, FileHandle, str): prmtop file, inpcrd file,
                                        tleap log text
        or:
            str: the tleap script (if script_only=True)

    '''
    _check_available('tleap')
    inputs = ['script', 'system.pdb']
    outputs = ['system.prmtop', 'system.inpcrd', 'STDOUT']
    script = "".join([f'source leaprc.{f}\n' for f in ff])

    if solvate:
        if solvate not in ['oct', 'box', 'cube']:
            raise ValueError(f'Error: unrecognised solvate option "{solvate}"')
        water_box = 'TIP3PBOX'
        for f in ff:
            if 'water.opc' in f:
                water_box = 'OPCBOX'
            elif 'water.tip4pew' in f:
                water_box = 'TIP4PEWBOX'
            elif 'water.tip4p2005' in f:
                water_box = 'TIP4P2005BOX'
            elif 'water.spce' in f:
                water_box = 'SPCEBOX'
    if het_names:
        if len(het_names) > 0:
            for r in het_names:
                _check_exists(Path(het_dir) / f'{r}.frcmod')
                _check_exists(Path(het_dir) / f'{r}.mol2')
                script += f'loadamberparams {r}.frcmod\n'
                script += f'{r} = loadmol2 {r}.mol2\n'
                inputs += [f'{r}.mol2', f'{r}.frcmod']

    if script_only:
        script += f"system = loadpdb {amberpdb}\n"
    else:
        script += "system = loadpdb system.pdb\n"
    if solvate == "oct":
        script += f"solvateoct system {water_box} {padding}\n"
    elif solvate == "cube":
        script += f"solvatebox system {water_box} {padding} iso\n"
    elif solvate == "box":
        script += f"solvatebox system {water_box} {padding}\n"
    if solvate is not None:
        script += f"addions system Na+ {n_na}\naddions system Cl- {n_cl}\n"

    script += "saveamberparm system system.prmtop system.inpcrd\nquit"
    if script_only:
        return script

    tleap = SubprocessTask('tleap -f script')
    tleap.set_inputs(inputs)
    tleap.set_outputs(outputs)
    fh = FileHandler()
    scriptfile = fh.create('scriptfile')
    scriptfile.write_text(script)
    args = [scriptfile, amberpdb]
    if het_names:
        if len(het_names) > 0:
            for r in het_names:
                args += [f'{Path(het_dir)}/{r}.mol2',
                         f'{Path(het_dir)}/{r}.frcmod']
    prmtop, inpcrd, stdout = tleap(*args)
    if prmtop is None or inpcrd is None:
        raise RuntimeError(f'Error in leap: {stdout}')
    return prmtop, inpcrd, stdout


def make_leap(inpdb, outinpcrd, outprmtop, het_names=None,
              het_dir='.',
              forcefields=None,
              solvate=None,
              ion_molarity=None,
              padding=10.0):
    """
    Generate a tleap script to prepare a system for MD simulation.

    Args:
        inpdb (path-like): name of input PDB file
        outinpcrd (path-like): name of output inpcrd file
        outprmtop (path-like): name of output prmtop file
        het_names (None or list): 3-letter residue names for heterogens
        forcefields (None or list): List of forcefields to use
        solvate (None or str): Solvation option - can be 'box',
                               'cube', or 'oct'.
        padding (float): minimum distance from any solute atom
                        to a periodic box boundary (Angstroms)
        ion_molarity (float or None): If not None, add Na+ and Cl- ions
                                      to reach this molarity (M)
        het_dir (str): Directory to search for heterogen parameter files
    Returns:
        str: The tleap script as a string
    """
    _check_exists(inpdb)
    if not forcefields:
        print('Warning: no forcefields specified, '
              'defaulting to "protein.ff14SB"', file=sys.stderr)
        forcefields = ['protein.ff14SB']
    has_protein_ff = False
    for ff in forcefields:
        if 'protein' in ff:
            has_protein_ff = True
    if not has_protein_ff:
        if 'gaff2' in forcefields:
            print('Warning: no protein forcefield specified,'
                  ' defaulting to protein.ff19SB.', file=sys.stderr)
            forcefields.append('protein.ff19SB')
        else:
            print('Warning: no protein forcefield specified,'
                  ' defaulting to protein.ff14SB.', file=sys.stderr)
            forcefields.append('protein.ff14SB')
    if solvate:
        if solvate not in ['oct', 'box', 'cube']:
            raise ValueError(f'Error: unrecognised solvate option "{solvate}"')
        water_ff = False
        for ff in forcefields:
            if 'water' in ff:
                water_ff = True
        if not water_ff:
            print('Warning: no water forcefield specified but'
                  ' solvation required.', file=sys.stderr)
            if 'protein.ff19SB' in forcefields:
                print('Defaulting to "water.opc" forcefield.',
                      file=sys.stderr)
                forcefields.append('water.opc')
            else:
                print('Defaulting to "water.tip3p" forcefield.',
                      file=sys.stderr)
                forcefields.append('water.tip3p')

    if het_names is not None:
        if 'gaff' not in forcefields and 'gaff2' not in forcefields:
            print('Warning - heterogens are present but no gaff/gaff2 '
                  'forcefield has been specified.', file=sys.stderr)
            print('Will default to using "gaff".', file=sys.stderr)
            forcefields.append('gaff')

    if not ion_molarity:
        try:
            script = leap(
                inpdb, forcefields, het_names=het_names,
                solvate=solvate, padding=padding, het_dir=het_dir,
                script_only=True)
            script = script.replace('system.prmtop', str(outprmtop))
            script = script.replace('system.inpcrd', str(outinpcrd))
            return script
        except RuntimeError as e:
            print(f'Error in leap:\n{e}')
            exit(1)

    try:
        prmtop, inpcrd, stdout = leap(
            inpdb, forcefields, het_names=het_names,
            solvate=solvate, padding=padding, het_dir=het_dir)
    except RuntimeError as e:
        print(f'Error in leap:\n{e}')
        exit(1)
    if ion_molarity:
        ttmp = mdt.load(inpcrd, top=prmtop)
        n_waters = len(ttmp.topology.select('name O and resname HOH'))
        n_na = len(ttmp.topology.select('name "Na+"'))
        n_cl = len(ttmp.topology.select('name "Cl-"'))
        n_ions = int(ion_molarity * n_waters/55.56)
        if n_na > 0:
            n_cl = n_ions - n_na
            n_na = n_ions
        else:
            n_na = n_ions - n_cl
            n_cl = n_ions
        n_na = max(n_na, 0)
        n_cl = max(n_cl, 0)

        try:
            script = leap(
                inpdb, forcefields, het_names=het_names,
                solvate=solvate, padding=padding, het_dir=het_dir,
                n_na=n_na, n_cl=n_cl, script_only=True)
            script = script.replace('system.prmtop', str(outprmtop))
            script = script.replace('system.inpcrd', str(outinpcrd))
            return script
        except RuntimeError as e:
            print(f'Error in leap:\n{e}')
            exit(1)


parser = ArgumentParser(description="Generate tleap input script from PDB")
parser.add_argument('--inpdb', help='Input PDB file', required=True)
parser.add_argument('--outinpcrd', help='Output AMBER .inpcrd file',
                    required=True)
parser.add_argument('--outprmtop', help='Output AMBER .prmtop file',
                    required=True)
parser.add_argument('--forcefields', nargs='*', help='Force fields to use')
parser.add_argument('--het_names', nargs='*',
                    help='Names of heterogen residues')
parser.add_argument('--solvate', help='Type of water box to use',
                    choices=['box', 'cube', 'oct'])
parser.add_argument('--padding',
                    help='minimum distance of solute atoms from box edge',
                    type=float, default=10.0)
parser.add_argument('--het_dir', help='Directory for heterogen files',
                    default='.')
parser.add_argument('--ion_molarity', type=float,
                    help='Target ionic strength (M)')

parser.add_argument('--version', action='version', version=__version__)

parsed_args = parser.parse_args()

try:
    result = make_leap(**vars(parsed_args))
except Exception as e:
    print(f'Error: {e}', file=sys.stderr)
    sys.exit(1)
print(result)
