## reduced_alphabet

Description of reduced_alphabet

### Compiling

Getting this to compile;

```

git clone https://github.com/seqan/seqan.git
git clone https://github.com/martinjvickers/reduced_alphabet.git
cd reduced_alphabet
cmake ../reduced_alphabet \
   -DCMAKE_MODULE_PATH=../seqan/util/cmake \
   -DSEQAN_INCLUDE_PATH=../seqan/include \
   -DCMAKE_CXX_FLAGS=-std=c++14 \
   -DCMAKE_BUILD_TYPE=Release \
   -DCMAKE_CXX_ARGS="-DSEQAN_DISABLE_VERSION_CHECK=YES" 
make
```

### Running

```
./reduced_alphabet -q example_data/protein.fasta -r example_data/yeast.fasta -o meh.text -k 3
```
