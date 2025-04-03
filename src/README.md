## Running

```bash
cargo run -- --query nl43 --distance 8 --assay random-reverse --file1 ./sim_data/mockseq_r1.fasta --file2 ./sim_data/mockseq_r2.fasta
```

```bash
./virust-splicing -q nl43 -d 8 -a random-reverse -1 ./sim_data/mockseq_r1.fasta -2 ./sim_data/mockseq_r2.fasta
```

<small>Note that inputs use '-' instead of '\_' , use assay "random-reverse" instead of "random_reverse"</small>
