# For fremvisning i videopresentasjon

> Laget for fremvisning i forbindelse med videopresentasjonen

## Triplet alignment

### Test-streng fra oppgaven i O2

>TTATGTTTTAAGGATGGGGCGTTAGTT
>
>--> TTT,GGGCGT

TTA TGT TTT AAG GAT GGG GCG TTA GTT
 TT **ATG** TTT *TAA* GGA TGG GGC GTT AGT T
 TT **ATG** TTT *TAA* GG **ATG** GGG CGT *TAG* TT

## Find genes algoritmen

TTATGTTTTAAGGATGGGGCGTTAGTT

1. Finn første start kodon
   >TT*ATG*TTTTAAGGATGGGGCGTTAGTT
   > `seq_start_pos: 2`
2. Iterer over triplene etter start kodon
   - Mellomlagre hver triplet
   - Bryt løkken hvis start triplet
   - Lagre gyldig gen og bryt løkken hvis
     stop triplet
   >TT **ATG** TTT *TAA* GGATGGGGCGTTAGTT

   ```python
   triplet: 'TAA'
   cur_seq: ['TTT']
   # Lagre gyldig gen
   valid_genes: ['TTT']
   ```

3. Finn neste start kodon
   >TT**ATG**TTT*TAA*GG **ATG** GGGCGTTAGTT

4. Hvis strict align og alignment er ugyldign
   >TT**ATG**TTT*TAA*GG **ATG** GGGCGTTAGTT

   ```python
   first_seq_start_pos: 2
   seq_start_pos: 13
   (13 - 2) % 3 = 2  # Ugyldig, 2 != 0
   ```

   - Finn neste start kodon
   - Gjenta 4.
5. Gjenta 2-4 til alle start kodon er funnet.
