<table align="center">
<tr><td align="center" width="10000">

<p>
    <img src = "./bullviso/assets/images/bv_network_diagram.png" width = "200">
</p>

# <strong> B U L L V I S O </strong>

<p>
    <a href="https://linkedin.com/in/conorrankine" > Dr. Conor D. Rankine </a> + <a href="https://york.ac.uk/chemistry/people/pmcgonigal/" > Prof. Paul McGonigal </a> <br> @ <a href="https://york.ac.uk" >The University of York </a>
</p>

<p>
    <a href="https://mcgonigalgroup.com/" > The McGonigal Group </a>
    <br>
    with special thanks to
    <br>
    <a href="https://linkedin.com/in/mariia-kuznetsova-9870bb217/" > Maria Kuznetsova </a>
    <br>
    <a href="https://linkedin.com/in/robives5" > Rob Ives </a>
    <br>
    <a href="https://linkedin.com/in/will-maturi-930738196" > Will Maturi </a>
</p>

<p>
If you enjoy BULLVISO, don't forget to cite:
</p>

<p>
    <a href="https://doi.org/10.1039/D4SC03700F" > A Guide to Bullvalene Stereodynamics </a>
    <br>
    <i>Chem. Sci.</i>, 2024, <b>15</b>, 14608-14617 (DOI: 10.1039/D4SC03700F).
</p>

<p>
    <a href="https://doi.org/10.1039/D4SC03699A" > Correlated Shapeshifting and Configurational Isomerization </a>
    <br>
    <i>Chem. Sci.</i>, 2024, <b>15</b>, 14618-14624 (DOI: 10.1039/D4SC03699A).
</p>

</td></tr></table>

#

## ⚙️ SETUP

Our favourite way to use BULLVISO is with <a href="https://docs.astral.sh/uv/" > uv </a> from <a href="https://astral.sh/" > Astral</a>! It's super-easy to install from your terminal (if you don't have it already) with:

```
pip install uv
```

You can use `uv` to create and activate a virtual environment:

```
uv venv bullviso --python 3.10
source bullviso/bin/activate
```

and then install BULLVISO directly into your new virtual environment with:

```
uv pip install git+https://www.github.com/conorrankine/bullviso
```

Don't forget to keep BULLVISO up to date; it's still under development, after all! It's easy to do with:

```
uv pip install --upgrade git+https://www.github.com/conorrankine/bullviso
```

When you're done with BULLVISO, you can deactivate your virtual environment with:

```
deactivate
```

and, on returning to the same directory where you created your virtual environment, reactivate it later with:

```
source bullviso/bin/activate
```

Now you're good to go!

## 🦾 QUICKSTART

### BASIC USAGE

BULLVISO has one core command line script: `bullviso`.

#### HOMOSUBSTITUTED BULLVALENES

##### MINIMUM-ENERGY GEOMETRIES

You can use `bullviso` to generate the unique constitutional isomers of a monosubstituted bullvalene by passing the SMILES string of the substituent as a command line argument, *e.g.*,

```
bullviso "C"
```

generates the four unique constitutional isomers of methylbullvalene.

The Cartesian coordinates for the minimum-energy geometries are written to `<output>/minima/`. Instructions on how to set the `<output>` directory and on the supported file formats are available under the section on [output options](#output-options).

You can use the `-n [N_SUBS]` flag to generate di- , tri-, *etc.* substituted (*i.e* *n*-substituted) bullvalenes, *e.g.*,

```
bullviso "C" -n 3
```

generates the 42 unique constitutional isomers of trimethylbullvalene.

For larger substituents, you might also want to specify an alternative attachment point on the SMILES string of the substituent. BULLVISO attachs the substituent *via* the first atom in the SMILES string by default.

You can use the `-a [SUB_ATTACH_IDX]` flag to specify an alternative attachment point, *e.g.*,

```
bullviso "CCCC"
```

generates the four unique constitutional isomers of *n*-butylbullvalene, while

```
bullviso "CCCC" -a 2
```

generates the four unique constitutional isomers of 2-butylbullvalene.

BULLVISO can accept SMILES strings for charged and radical substituents, *e.g.*,

```
bullviso "[O-]"
```

generates the four unique constitutional isomers of the bullvalene alcoholate anion (the deprotonated form of hydroxybullvalene), while

```
bullviso "[NH3+]"
```

generates the four unique constitutional isomers of the bullvalene ammonium cation (the protonated form of aminobullvalene).

BULLVISO can also accept SMILES strings with explicit stereochemical specifications for atoms (*i.e.* *R*/*S* chirality) and bonds (*i.e.* *E*/*Z* stereochemistry), *e.g.*,

```
bullviso "N[C@@H](C)C(=O)" -a 4
```

and

```
bullviso "N[C@H](C)C(=O)" -a 4
```

output the four unique constitutional isomers of (*L*-alanyl)bullvalene and (*D*-alanyl)bullvalene, respectively.

⚠️ **Note**: the H, CH<sub>3</sub>, and C(=O)R groups appear clockwise (@@) and counterclockwise (@) when viewed down the N–C bond for (*L*-alanyl)bullvalene and (*D*-alanyl)bullvalene, respectively. The attachment point, set using the `-a [SUB_ATTACH_IDX]` flag, is to the **fourth** (**not the fifth!**) atom of the substituent since explicit hydrogens in SMILES strings are skipped, even in cases where they inform the stereochemistry.

Similarly,

```
bullviso "C(=C/F)\F"
```

and

```
bullviso "C(=C/F)/F"
```

output the four unique constitutional isomers of (*E*-1,2-difluoroethenyl)bullvalene and (*Z*-1,2-difluoroethenyl)bullvalene, respectively.

##### TRANSITION STATE GEOMETRIES

In addition to the unique minimum-energy geometries, BULLVISO also generates the transition-state (TS) geometries that connect the minima *via* Cope rearrangement.

The Cartesian coordinates for the TS geometries are written to `<output>/transition_states/`.

If you only need the minimum-energy geometries, you can disable TS geometry generation with the `-no-ts` flag, *e.g.*,

```
bullviso "C" -no-ts
```

generates the four unique constitutional isomers of methylbullvalene *without* the TS geometries that connect the minima *via* Cope rearrangement.

#### CONFORMER GENERATION

BULLVISO has an embedded multi-step sequential workflow for the generation, geometry optimisation, and selection of a diverse set of low-energy molecular conformations.

Conformers are:

1. **embedded** using the ETKDGv3 distance geometry algorithm;

2. **optimised** using a calculator implementing either the Merck Molecular Forcefield (MMFF), Universal Forcefield (UFF), or extended tight binding (XTB) semiempirical quantum mechanical framework;

3. **filtered** to remove high-energy instances that exceed a cutoff energy threshold relative to the lowest-energy conformer;

4. **clustered** using the Butina algorithm to group structurally similar instances based on pairwise root-mean-squared distance (RMSD);

5. **selected** from the Butina clusters, retaining only the lowest-energy instance belonging to each;

6. **sorted** in ascending order by energy.

##### CONTROLLING CONFORMER GENERATION

The conformation generation workflow runs '*under the hood*', and BULLVISO outputs only the lowest-energy conformational isomer for each unique constitutional isomer by default but, if you're working with large or floppy/flexible substituents, you might want to output many conformational isomers for each unique constitutional isomer (ahead of, *e.g.*, population analysis).

You can use the `-m [M_CONFS]` flag to output the *m* lowest-energy conformational isomers of each unique constitutional isomer instead, *e.g.*,

```
bullviso "CCCC" -n 2 -m 6
```

outputs (up to a maximum of) the six lowest-energy conformational isomers for each of the 15 unique constitutional isomers of dibutylbullvalene.

The workflow is '*black-box*' and you don't *have* to alter the presets, although it's easy to do this from the command line if you need to configure/customise it! The relevant command line arguments are summarised below:

| Argument                 | Flag  | Description                                                             | Default |
|:-------------------------|:------|:------------------------------------------------------------------------|--------:|
| `--embed_n_confs`        | `-en` | maximum number of conformational isomers to embed                       | None    |
| `--embed_rmsd_threshold` | `-er` | RMSD threshold (Angstroem) for deduplicating embeddings                 | 0.1     |
| `--embed_timeout`        | `-et` | timeout (seconds) for conformer embedding                               | None    |
| `--embed_seed`           | `-es` | random seed for conformer embedding                                     | None    |
| `--calculator_type`      | `-c`  | calculator type for conformer optimisation                              | 'MMFF'  |
| `--max_iter`             | `-it` | maximum number of iterations for conformer optimisation                 | 600     |
| `--energy_threshold`     | `-e`  | energy threshold (kcal/mol) for energy-based conformer filtration       | 10.0    |
| `--rmsd_threshold`       | `-r`  | RMSD threshold (Angstroem) for distance-based conformer clustering      | 0.5     |
| `--n_proc`               | `-np` | number of (parallel) processes for conformer embedding and optimisation | 1       |

##### XTB SUPPORT

BULLVISO supports conformer optimisation and energy calculation at the XTB level through an interface to the <a href=https://github.com/grimme-lab/xtb > `xtb` </a> binary from the <a href=https://www.chemie.uni-bonn.de/grimme > Grimme Group </a>). To use XTB in BULLVISO, you'll need to install the standalone XTB program and ensure that the `xtb` binary is accessible through your system `$PATH`.

You can then use the `-c [CALCULATOR_TYPE]` flag to select the XTB calculator, *e.g.*,

```
bullviso "c1ccc(F)cc1" -n 2 -m 10 -c "xtb" 
```

outputs (up to a maximum of) the ten lowest-energy conformational isomers for each of the 15 unique constitutional isomers of di(*p*-fluorophenyl)bullvalene, optimised, filtered, and sorted by energy at the XTB level.

BULLVISO uses the <a href=https://doi.org/10.1021/acs.jctc.8b01176 > GFN2-xTB </a> parameterisation scheme by default. It is not currently possible to select an alternative parameterisation scheme.

⚠️ **Note**: conformer optimisation and energy calculation at the XTB level is more time- and resource-intensive than at the MMFF/UFF level but provides considerably better geometries and energies that often approximate the results of higher-level electronic structure calculations [*e.g.* density functional theory (DFT)] well. If you need reliable geometries and energies when working with, *e.g.*, electron-rich, exotic, or interacting (*e.g.*, hydrogen bonding) substituents, XTB is your friend!

#### HETEROSUBSTITUTED BULLVALENES

You can generate heterosubstituted bullvalenes (*i.e.* bullvalenes with more than one kind of substituent) by passing the SMILES strings of multiple different substituents as separate command line arguments, *e.g.*,

```
bullviso "C" "CC" "CCC"
```

generates the 240 unique constitutional isomers of (methyl,ethyl,*n*-propyl)bullvalene.

You can use `-n [N_SUBS]` and `-a [SUB_ATTACH_IDX]` flags to modify the number of substituents and/or the attachment points on the SMILES strings of the substituents by passing lists of values as command line arguments where each element of the list corresponds to each of the SMILES strings defining the substituents, *e.g.*,

```
bullviso "C" "CC" "CCC" -n [2,2,1]
```

generates the 2520 unique constitutional isomers of (dimethyl,diethyl,*n*-propyl)bullvalene, while

```
bullviso "C" "CC" "CCC" -n [2,2,1] -a [1,1,2]
```

generates the 2520 unique constitutional isomers of (dimethyl,diethyl,isopropyl)bullvalene.

#### OUTPUT OPTIONS

The parent output directory (`<output>`) contains two subdirectories, `<output>/minima` and `<output>/transition_states`, containing the minimum-energy and TS geometries, respectively. Beneath these subdirectories, BULLVISO organises the outputted geometries into separate subdirectories labelled using the bullvalene isomer barcode and, beneath these, into further subdirectories labelled using the bullvalene isomer barcode and conformational isomer number.

A simple schema for the output directory hierarchy is given here:

```
<output>/
├── minima/
│   ├── 0000010010/
│   │   ├── 0000010010_001/
│   │   │   └── 0000010010_001.xyz
│   │   ├── 0000010010_002/
│   │   │   └── 0000010010_002.xyz
│   │   └── …
│   └── …
└── transition_states/
    ├── 0000010010/
    │   └── …
    └── …
```


By default, the parent output directory (`<output>`) is created in the present working directory (PWD), although you can use the `-o [OUTPUT_DIR]` flag to create the output directory hierarchy somewhere else on your system by providing an alternative path. If the output directory doesn't exist, BULLVISO will try to create it.

BULLVISO outputs the structures of the substituted bullvalenes in .xyz format by default, although you can use the `-f [OUTPUT_FILETYPE]` flag to output in an alternative format, *e.g.* .sdf/.mol (`-f sdf`), or a basic customisable input file for <a href="https://gaussian.com/gaussian16/" > Gaussian </a> (`-f gaussian`) or <a href="https://kofo.mpg.de/en/research/services/orca" > Orca </a> (`-f orca`).

### HELP

Don't forget that you can always get a reminder of, and help with, the available command line flags by entering:

```
bullviso --help
```

### ADVANCED USAGE:

#### BRIDGING SUBSTITUENTS

#### DOUBLY-BRIDGING SUBSTITUENTS

It's also possible to use the `-a [SUB_ATTACH_IDX]` flag to specify multiple attachment points on the same substituent and create a '*bridging*', or '*bifunctional*', group *e.g.*,

```
bullviso "CC(N)=O" -a [[1,3]]
```

generates the 30 unique constitutional isomers of (lactam)bullvalene. BULLVISO won't output all 30 in this case, though - some of these unique constitutional isomers are mechanically impossible (*e.g.* isomers where the lactam bridges between the $\alpha$ and $\delta$ positions of the bullvalene).

Note the double set of square brackets around the argument passed with the `-a [SUB_ATTACH_IDX]` flag (`[[1,3]]`) that turn it into a nested (sub)list: this is to specify that the first substituent (acetamide; SMILES = "CC(N)=O"), is attached to the bullvalene at the atomic indices 1 and 3. A single set of square brackets around the argument passed with the `-a [SUB_ATTACH_IDX]` flag (`[1,3]`) would specify a first substituent attached at atomic index 1 and a second substituent attached at atomic index 3 with reference to their respective SMILES strings; BULLVISO would throw an error here as only a single SMILES string is specified.

You can generate substituted bullvalenes with both singly- and doubly-attached substituents using the same syntax, *e.g.*,

```
bullviso "C" "CC(N)=O" -n [2,1] -a [1,[1,3]]
```

generates the 840 unique constitutional isomers of (dimethyl,lactam)bullvalene.

#### MULTIPLY-BRIDGING/MULTIPLEXING SUBSTITUENTS

As long as the total number of attachment points remains below nine, there are no limits as to how many attachment points can be specified for a single substituent, allowing you to create '*multiply-bridging*', '*multiplexed*', or '*polyfunctional*' groups - as long as your ideas are mechanically sound!

#### SYMMETRIC VS. ASYMMETRIC BRIDGING SUBSTITUENTS

BULLVISO considers each attachment point on a single substituent as unique, *i.e.*, for a substituted bullvalene with a single substituent attached at two points, like (lactam)bullvalene, the [0000000012] and [0000000021] constitutional isomers are considered inequivalent.

This is fine for (lactam)bullvalene since the [1] and [2] bits in the bullvalene isomer barcode correspond to the inequivalent attachment points at the atomic indices 1 (C) and 3 (N), respectively, on the substituent (acetamide; SMILES = "CC(N)=O"). For symmetric substituents attached at two (or more) points, however, this assumption is *not* necessarily valid, *e.g.*,

```
bullviso "CSSC" -a [[1,4]]
```

generates 30 constitutional isomers for (disulfide)bullvalene, but only 15 of these are unique since the [1] and [2] bits in the bullvalene isomer barcode are interchangable (accounting for the symmetry of the substituent): *e.g.*, the structures of the [0000000012] and [0000000021] constitutional isomers are identical.

## 🗒️ LICENSE

BULLVISO is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License (GPL) as published by the Free Software Foundation, either Version 3 of the license (GPLv3), or (at your option) any later version.

BULLVISO is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the Lesser GNU GPL for more details.

Unless you explicitly state otherwise, any contribution you intentionally submit for inclusion in BULLVISO is licensed as in the Lesser GNU GPL without any additional terms or conditions.
