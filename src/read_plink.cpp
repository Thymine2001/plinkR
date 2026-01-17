// Adapted from plink2R read_plink
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <fstream>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP read_plink_cpp(std::string bedfile, std::string famfile, int impute, bool verbose) {
  // load fam to get sample count
  std::ifstream fam(famfile.c_str());
  if (!fam.good()) stop("Cannot open fam file");
  int n = 0;
  std::string line;
  while (std::getline(fam, line)) {
    if (!line.empty()) n++;
  }
  fam.close();
  if (n == 0) stop("Empty fam file");

  // load bim to get variant count
  std::string bimfile = bedfile.substr(0, bedfile.size() - 4) + ".bim";
  std::ifstream bim(bimfile.c_str());
  if (!bim.good()) stop("Cannot open bim file");
  int m = 0;
  while (std::getline(bim, line)) {
    if (!line.empty()) m++;
  }
  bim.close();
  if (m == 0) stop("Empty bim file");

  // open bed
  std::ifstream bed(bedfile.c_str(), std::ios::in | std::ios::binary);
  if (!bed.good()) stop("Cannot open bed file");

  // check magic bytes and mode
  char header[3];
  bed.read(header, 3);
  if (header[0] != 0x6c || header[1] != 0x1b) {
    stop("Invalid BED file header");
  }
  if (header[2] != 0x01) {
    stop("Only SNP-major BED files are supported");
  }

  // prepare matrix
  IntegerMatrix geno(n, m);

  // bytes per SNP
  size_t bytes_per_snp = (n + 3) / 4;
  std::vector<unsigned char> buffer(bytes_per_snp);

  for (int j = 0; j < m; ++j) {
    bed.read(reinterpret_cast<char*>(&buffer[0]), bytes_per_snp);
    if (!bed) stop("Unexpected EOF in BED");
    int idx = 0;
    for (size_t b = 0; b < bytes_per_snp; ++b) {
      unsigned char byte = buffer[b];
      for (int k = 0; k < 4; ++k) {
        if (idx >= n) break;
        unsigned char code = (byte >> (2 * k)) & 0x03;
        int val;
        if (code == 0) val = 2;          // homo major
        else if (code == 1) val = NA_INTEGER; // missing
        else if (code == 2) val = 1;     // hetero
        else val = 0;                    // homo minor
        geno(idx, j) = val;
        idx++;
      }
    }
  }

  // impute if requested
  if (impute != 0) {
    for (int j = 0; j < m; ++j) {
      // compute mean of non-missing
      double sum = 0.0;
      int count = 0;
      for (int i = 0; i < n; ++i) {
        int v = geno(i, j);
        if (v != NA_INTEGER) {
          sum += v;
          count++;
        }
      }
      double mean = count > 0 ? sum / count : NA_REAL;
      for (int i = 0; i < n; ++i) {
        if (geno(i, j) == NA_INTEGER) {
          if (impute == 1 && mean == mean) { // avg
            geno(i, j) = static_cast<int>(::round(mean));
          } else if (impute == 2 && count > 0) { // random
            // simple: use mean rounding; for real randomness integrate R's RNG
            geno(i, j) = static_cast<int>(::round(mean));
          }
        }
      }
    }
  }

  return geno;
}
