#include <stdio.h>
#include <string.h>
#include "pwlib.h"

/**
 * Given an alignment ::transcript returns the length of its opseq on the
 * "from" or "to" sequence.
 *
 * @param tx The transcript of interest.
 * @param on A single character of either 'S' or 'T'.
 * @return A non-negative integer corresponding to the length of the indicated
 *   sequence that is covered by the transcript and -1 if an error occurs.
 */
int tx_seq_len(transcript* tx, char on) {
  if (on != 'S' && on != 'T') {
    return -1;
  }
  int sum = 0;
  for (int i = 0; i < strlen(tx->opseq); i++) {
    if ((on == 'S' && tx->opseq[i] != 'I') ||
        (on == 'T' && tx->opseq[i] != 'D')) {
      sum ++;
    }
  }
  return sum;
}


/**
 * Helper method for debugging purposes. Only works on system with procfs.
 */
void _print_mem_usage() {
  int tot[1];
  int res[1];
  FILE *fp = fopen("/proc/self/statm", "r");
  fscanf(fp, "%d %d", tot, res);
  printf("%.2f MB (tot), %.2f MB (res)\n",
    tot[0]/1024.0, res[0]/1024.0);
  fclose(fp);
}
