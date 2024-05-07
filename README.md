<table align="center">
<tr><td align="center" width="10000">

# <strong> B U L L V I S O </strong>

<p>
    <a href="https://linkedin.com/in/conorrankine"> Dr. Conor D. Rankine </a> and <a href="https://linkedin.com/in/conorrankine"> Dr. Paul McGonigal </a> @ <a href="https://york.ac.uk">The University of York </a>
</p>

<p>
    <a href="https://www.mcgonigalgroup.com/"> The McGonigal Group </a>
</p>

</td></tr></table>

#

## SETUP

If you have a working Python installation on your machine, the quickest way to get started with BULLVISO is to pop open a terminal and clone this repository:

```
git clone https://www.gitlab.com/conorrankine/bullviso
```

...change directory:

```
cd ./bullviso
```

...and install BULLVISO with pip:

```
pip install .
```

Now you're good to go!

## QUICKSTART

BULLVISO has one core command line script: `bullviso`. You can use `bullviso` to generate the unique constitutional isomers of a monosubstituted bullvalene by passing the SMILES string of the substituent as a command line argument, *e.g.*,

```
bullviso C
```

generates the four unique constitutional isomers of methylbullvalene.

You can use the `-n [N_SUBS]` flag to generate di- , tri-, *etc.* substituted (*i.e* *n*-substituted) bullvalenes, *e.g.*,

```
bullviso C -n 3
```

generates the 42 unique constitutional isomers of trimethylbullvalene.

For larger substituents, you might want to generate the conformational and/or configurational isomers of each constitutional isomer. BULLVISO generates only the lowest-energy conformational/configurational isomer (minimised under the Universal Force Field) of each constitutional isomer by default.

You can use the `-m [M_CONFS]` flag to generate the *m* lowest-energy conformational/configurational isomers of each constitutional isomer instead, *e.g.*,

```
bullviso CCCC -n 2 -m 6
```

generates (up to a maximum of) the six lowest-energy conformational/configurational isomers of each of the 15 unique constitutional isomers of dibutylbullvalene.

For larger substituents, you might also want to specify an alternative attachment point on the SMILES string of the substituent. BULLVISO attachs the substituent *via* the first atom in the SMILES string by default.

You can use the `-a [SUB_ATTACH_IDX]` flag to specify an alternative attachment point, *e.g.*,

```
bullviso CCCC
```

generates the four unique constitutional isomers of *n*-butylbullvalene, while

```
bullviso CCCC -a 2
```

generates the four unique constitutional isomers of 2-butylbullvalene.

BULLVISO outputs the structures of the substituted bullvalenes in .xyz format by default, although you can use the `-o [OUT_F_TYPE]` flag to output in an alternative format, *e.g.* a basic (customisable) input file for Gaussian (`-o gaussian`) or Orca (`-o orca`). 

Don't forget that you can always get a reminder of, and help with, the available command line flags by entering:

```
bullviso -h
```

## LICENSE

BULLVISO is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License (GPL) as published by the Free Software Foundation, either Version 3 of the license (GPLv3), or (at your option) any later version.

BULLVISO is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the Lesser GNU GPL for more details.

Unless you explicitly state otherwise, any contribution you intentionally submit for inclusion in BULLVISO is licensed as in the Lesser GNU GPL without any additional terms or conditions.