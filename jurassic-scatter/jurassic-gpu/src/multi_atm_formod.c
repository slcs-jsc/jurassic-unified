#include <string.h> // memcpy
#include <omp.h>
#include "jurassic.h" // ctl_t, obs_t, atm_t, read_*, write_*, formod_*PU

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

const int MAX_LIST_SIZE = 100;
const int MAX_FILENAME_LENGTH = 30;
const int MAX_OBS = 100;

int read_atm_list(char const *dirname, char const *filename, 
                   char atm_list[MAX_LIST_SIZE][MAX_FILENAME_LENGTH]) {
	FILE *in;
	char line[LEN], *tok, *saveptr;
  int i = 0;
	printf("Read atm list: %s/%s\n", dirname, filename);
	// Open file
	in = jur_mkFile(dirname, filename, "r");
	// Read line
	while (fgets(line, LEN, in)) {
		// Read data
		TOK_FIVE_ARGS(line, tok, "%s", atm_list[i++], &saveptr);
	}
	// Close file
	fclose(in);
	// Check number of points
	if(i < 1) ERRMSG("Could not read any data!");
	printf("Read atm list, found %d files\n", i);
  return i;
}

int read_atm_id(char const *dirname, char const *filename, 
                 int atm_id[MAX_OBS]) {
	FILE *in;
	char line[LEN], *tok, *saveptr;
  int i = 0;
	printf("Read atm list: %s/%s\n", dirname, filename);
	// Open file
	in = jur_mkFile(dirname, filename, "r");
	// Read line
	while (fgets(line, LEN, in)) {
		// Read data
		TOK_FIVE_ARGS(line, tok, "%d", atm_id[i++], &saveptr);
	}
	// Close file
	fclose(in);
	// Check number of points
	if(i < 1) ERRMSG("Could not read any data!");
	printf("Read atm_id, found %d indices\n", i);
  return i;
}


int main(
	 int argc,
	 char *argv[]) {
  
  static ctl_t ctl; // why static?

  atm_t *atm;

  obs_t *obs;
  
  /* Check arguments... */
  if (argc < 6)
    ERRMSG("Give parameters: <ctl> <obs> <atm_list> <atm_id> <rad>");

  ALLOC(obs, obs_t, 1);
  
  /* Read control parameters... */
  jur_read_ctl(argc, argv, &ctl);
  
  /* Read observation geometry... */
  jur_read_obs(".", argv[2], &ctl, obs);

  /* Read atmospheric data... */
  // jur_read_atm(".", argv[3], &ctl, atm);
  // to verify that we read in the atms correctly
  //   write_atm(".", "atm_ref.tab", &ctl, atm);
  //
  char atm_list[MAX_LIST_SIZE][MAX_FILENAME_LENGTH];
  int num_of_atms = read_atm_list(".", argv[3], atm_list);

  int atm_id[MAX_OBS];
  int atm_id_list_length = read_atm_id(".", argv[3], atm_id);

  /* Allocate... */
  ALLOC(atm, atm_t, num_of_atms);
  for(int i = 0; i < num_of_atms; i++)
    jur_read_atm(".", atm_list[i], &ctl, &atm[i]);


  /* Call forward model... */
  jur_formod(&ctl, &atm[0], obs);

  /* Write radiance data... for reference */
  jur_write_obs(".", argv[5], &ctl, obs);
  
  free(atm);
  free(obs);
  
  return EXIT_SUCCESS;
}
