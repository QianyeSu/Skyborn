#ifndef CDO_OPTIONS_H
#define CDO_OPTIONS_H

#include <vector>
#include <string>
#include "mpim_grid/grid_options.h"

#ifdef HAVE_CONFIG_H
#include "config.h" /* _FILE_OFFSET_BITS influence off_t */
#endif
#ifndef VERSION
#define VERSION "0.0.1"
#endif

namespace cdo
{
extern const char *progname;
extern const char *Version;
extern std::string FileSuffix;
extern bool stdinIsTerminal;
extern bool stdoutIsTerminal;
extern bool stderrIsTerminal;
}  // namespace cdo

enum struct MemType
{
  Native,
  Float,
  Double
};

namespace Options
{
extern long coresize;
extern int numStreamWorker;
extern int nsb;  // Number of significant bits
extern bool benchmark;
extern bool silentMode;
extern bool test;
extern bool fast;
extern bool lazyGridLoad;
extern bool force;

extern std::string filterSpec;

extern int cdoShuffle;
extern bool cdoCompress;
extern int cdoCompType;
extern int cdoCompLevel;
extern bool cdoInteractive;
extern bool cdoVerbose;
extern int CDO_Rusage;
extern bool cdoProcessInfo;
extern int cdoExitStatus;
extern bool Timer;

extern bool CheckDatarange;

extern int CDO_flt_digits;
extern int CDO_dbl_digits;

extern bool Use_FFTW;
extern bool VersionInfo;
extern int CMOR_Mode;

extern bool CDO_diagnostic;

extern MemType CDO_Memtype;

extern bool CDO_Async_Read;
extern bool CDO_task;

extern int CDO_Reduce_Dim;
extern int CDO_Append_History;
extern bool CDO_Reset_History;
extern bool PrintFilename;

extern unsigned Random_Seed;

extern int cdoChunkType;
extern int cdoChunkSize;
extern int cdoChunkSizeDimT;
extern int cdoChunkSizeDimZ;
extern int cdoChunkSizeDimY;
extern int cdoChunkSizeDimX;
extern bool cdoOverwriteMode;
extern bool cdoParIO;
extern bool cdoRegulargrid;
extern std::string cdoQueryParameter;
extern std::vector<std::string> cdoVarnames;
size_t cdo_num_varnames();

extern bool RemapGenerateWeights;

}  // namespace Options

namespace Threading
{
extern int ompNumMaxThreads;
extern int ompNumUserRequestedThreads;
extern bool cdoLockIO;
}  // namespace Threading

const char *cdo_comment(void);

void set_compression(int fileID, int filetype);

void cdo_set_search_radius(double searchRadius);
double cdo_get_search_radius(void);
void cdo_print_attributes(std::FILE *fp, int cdiID, int varID, int nblanks);

#endif
