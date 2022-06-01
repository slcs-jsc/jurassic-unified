#include <string.h> // memcpy
#include <omp.h>
#include "jurassic.h" // ctl_t, obs_t, atm_t, read_*, write_*, formod_*PU

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

const int MAX_LIST_SIZE = 100;
const int MAX_FILENAME_LENGTH = 30;
const int MAX_OBS = 100;

int read_name_list(char const *dirname, char const *filename, 
    char name_list[MAX_LIST_SIZE][MAX_FILENAME_LENGTH]) {
  FILE *in;
  char line[LENMAX], *tok, *saveptr;
  int i = 0;
  printf("Read name list: %s/%s\n", dirname, filename);
  // Open file
  in = jur_mkFile(dirname, filename, "r");
  // Read line
  while (fgets(line, LENMAX, in)) {
    // Read data
    TOK_FIVE_ARGS(line, tok, "%s", name_list[i++], &saveptr);
  }
  // Close file
  fclose(in);
  // Check number of points
  if(i < 1) ERRMSG("Could not read any data!");
  printf("Read name list, found %d files\n", i);
  return i;
}

int read_atm_id(char const *dirname, char const *filename, 
    int32_t atm_id[MAX_OBS]) {
  FILE *in;
  char line[LENMAX], *tok, *saveptr;
  int i = 0;
  printf("Read atm list: %s/%s\n", dirname, filename);
  // Open file
  in = jur_mkFile(dirname, filename, "r");
  // Read line
  while (fgets(line, LENMAX, in)) {
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

// at the moment we use only CPU version
void formod_CPU(ctl_t const *ctl, atm_t *atm, obs_t *obs,
    int32_t const *atm_id, aero_t const *aero); 

int main(
    int argc,
    char *argv[]) {

  static ctl_t ctl; // why static?

  atm_t *atm;

  obs_t *obs;

  char atm_list[MAX_LIST_SIZE][MAX_FILENAME_LENGTH];
  char obs_list[MAX_LIST_SIZE][MAX_FILENAME_LENGTH];
  char rad_list[MAX_LIST_SIZE][MAX_FILENAME_LENGTH];

  int32_t atm_id[MAX_OBS];

  if (argc < 6)
    ERRMSG("Give parameters: <ctl> <obs> <atm_list> <atm_id> <rad>");

  jur_read_ctl(argc, argv, &ctl);

  int num_of_atms = read_name_list(".", argv[3], atm_list);

  int atm_id_list_length = read_atm_id(".", argv[4], atm_id);

  int num_of_obs = read_name_list(".", argv[2], obs_list);

  int num_of_rads = read_name_list(".", argv[5], rad_list);

  assert(num_of_obs == num_of_rads);

  ALLOC(atm, atm_t, num_of_atms);
  for(int i = 0; i < num_of_atms; i++)
    jur_read_atm(".", atm_list[i], &ctl, &atm[i]);

  ALLOC(obs, obs_t, num_of_obs);
  for(int i = 0; i < num_of_obs; i++)
    jur_read_obs(".", obs_list[i], &ctl, &obs[i]);

  int total_num_of_rays = 0;
  for(int i = 0; i < num_of_obs; i++)
    total_num_of_rays += obs[i].nr;

  printf("list length vs. obs->nr: %d vs. %d\n", atm_id_list_length, total_num_of_rays);
  assert(atm_id_list_length == total_num_of_rays);
  for(int i = 0; i < total_num_of_rays; i++) 
    assert(atm_id[i] < num_of_atms);

  // FIXME: without this line there is a problem with barrier inside get_tbl(..) function
  // FIXME: also, SIGSEGV in useGPU case, but only when using gdb :o
  jur_table_initialization(&ctl); 

  jur_formod_multiple_packages(&ctl, atm, obs, num_of_obs, atm_id, NULL);

  for(int i = 0; i < num_of_rads; i++) {
    jur_write_obs(".", rad_list[i], &ctl, &obs[i]);
  }
  
  // ------------------------
  // testing muti - single - multi - single scenario

  atm_t *one_atm;
  ALLOC(one_atm, atm_t, 1);
	memcpy(one_atm, &atm[0], sizeof(atm_t));
  obs_t *one_obs;
  ALLOC(one_obs, obs_t, 1);
  memcpy(one_obs, &obs[0], sizeof(obs_t));

  jur_formod_multiple_packages(&ctl, one_atm, one_obs, 1, NULL, NULL);
  jur_formod_multiple_packages(&ctl, atm, obs, num_of_obs, atm_id, NULL);
  jur_formod_multiple_packages(&ctl, one_atm, one_obs, 1, NULL, NULL);
  jur_formod_multiple_packages(&ctl, atm, obs, num_of_obs, atm_id, NULL);

  free(atm);
  free(obs);

  return EXIT_SUCCESS;
}
