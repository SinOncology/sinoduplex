
<p align="center">
  <a href="http://www.sinotechgenomics.com">
    <img height="70" src="http://www.sinotechgenomics.com/Upload/0/WebsiteLogo/WebsiteLogo_20170620092534731.png">
  </a>
  <h1 align="center">SinoDuplex</h1>
</p>


SinoDuplex is an improved duplex consensus sequence generation module. It increased the output of duplex sequencing technology and thus made it more cost effective.

# Installing SinoDuplex
-----------
Download the source code of SinoDuplex.

`git clone https://github.com/SinOncology/sinoduplex.git`
`cd sinoduplex`

Edit the referenc-data.txt to define the folder contaning nib file and file containing barcodes
Then run the following shell script to compile SinoDuplex

`./install_sinoduplex.sh`

The binary of SinoDuplex is in sinoduplex/bin/sinotools 

# Running SinoDuplex
--------

`sinoduplex/bin/sinotools duplex -i <InputBamFile> -o <OutputPrefix> [-c <chrom>][-s] `

# Citation
-----

SinoDuplex: an improved duplex sequencing approach to detect low-frequency variants in plasma cfDNA samples. Yongzhe Ren et.al, submitted. 