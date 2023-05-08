/* Copyright Jérémie Burgalat (2010-2015,2017)
 * 
 * jeremie.burgalat@univ-reims.fr
 * 
 * This software is a computer program whose purpose is to provide configuration 
 * file and command line arguments parsing features to Fortran programs.
 * 
 * This software is governed by the CeCILL-B license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL-B
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 * 
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability. 
 * 
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or 
 * data to be ensured and,  more generally, to use and operate it in the 
 * same conditions as regards security. 
 * 
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-B license and that you accept its terms.
 */



/* ------------------ WRAPPERS METHODS FOR FORTRAN BINDINGS ------------------ */

/**
 * Gets the current umask (in decimal system)
 * @return An int with the current umask set for the current session
 */
int c_umask();

/** 
 * Get directory name of input path 
 * @param[in] in A C string with the input path
 * @return A pointer to a char array with the directory name of @bti{input} path.
 * @note On error, a NULL pointer is returned.
 * @warning In any case, the returned pointer must be freed in the Fortran counterpart 
 * (using fsystem::free_c).
 */
char * c_dirname(const char *in);

/** 
 * Get base name of input path 
 * @param[in] in A C string with the input path
 * @return A pointer to a char array with the base name of @bti{input} path.
 * @note On error, a NULL pointer is returned.
 * @warning In any case, the returned pointer must be freed in the Fortran counterpart 
 * (using fsystem::free_c).
 */
char* c_basename(const char *in);

/** 
 * Get the current working directory.
 * @return A pointer to a char array with the current workind directory.
 * @note On error, a NULL pointer is returned.
 * @warning In any case, the returned pointer must be freed in the Fortran counterpart 
 * (using fsystem::free_c).
 */
char* c_getcwd();


/** 
 * Get the realpath of input path.
 * @param[in] input A C string with the input path
 * @return A pointer to a char array with the realpath of @bti{input} path.
 * @note On error, a NULL pointer is returned.
 * @warning In any case, the returned pointer must be freed in the Fortran counterpart 
 * (using fsystem::free_c).
 */
char* c_realpath(const char *input);

/**
 * Get the relative path of two file paths
 * @details The method computes the relative path of can_fname to can_reldir if
 * possible.
 * @param path string with the path to compute in relative representation
 * @param reldir a directory path from which output should be relative to
 * @return A pointer to a char array with the relative path.
 * @note On error, a NULL pointer is returned.
 * @warning In any case, the returned pointer must be freed in the Fortran counterpart 
 * (using fsystem::free_c).
 */
char* c_relpath(const char *path, const char *reldir) ;

/**
 * Get the corresponding name of the given user id
 * @param[in] uid An integer with a user id
 * @return A pointer to a char array with the user name.
 * @note On error, a NULL pointer is returned.
 * @warning In any case, the returned pointer must be freed in the Fortran counterpart 
 * (using fsystem::free_c).
 */
char* c_uname(int uid);

/**
 * Get the corresponding name of the given group id
 * @param[in] gid An integer with a group id
 * @return A pointer to a char array with the group name.
 * @note On error, a NULL pointer is returned.
 * @warning In any case, the returned pointer must be freed in the Fortran counterpart 
 * (using fsystem::free_c).
 */
char* c_gname(int gid);

/**
 * Get last errno from C library
 * @return An int with last errno saved in C library.
 */
int c_get_errno();

/** 
 * Get the error message of the given error id
 * @param err An integer with the error id
 * @return A pointer to a char array with the group name.
 * @note On error, the hard-coded message "Unknown error" is returned.
 * @warning In any case, the returned pointer must be freed in the Fortran counterpart 
 * (using fsystem::free_c).
 */
char* c_strerror(int err);

/**
 * Creates directories recursively
 * @note The method is simple stripped copy of mkdir tool source code.
 * @param[in] path A C string with the path of the directory to create
 * @param[in] mode A mode_t structure with the permission of the directory
 * @return An integer with 0 on success, last errno on failure
 */
int c_mkdirp(const char *path, mode_t mode); 

/**
 * Rename a path
 * @param old A string with the (valid) path to rename
 * @param new A string with the new name of the path
 * @return An integer with 0 on success, last errno on failure
 */
int c_rename(const char *old, const char *new);

/**
 * Change path permissions
 * @param path A string with the path 
 * @param mode A integer with the new permissions to set
 * @return An integer with 0 on success, last errno on failure
 */
int c_chmod(const char *path, mode_t mode);

/**
 * Change current working directory
 * @param path A C string with the new path
 * @return An integer with 0 on success, last errno on failure
 */
int c_chdir(const char *path);

/**
 * Create directory
 * @param[in] path A C string with the path of the directory to create
 * @param[in] mode A mode_t structure with the permission of the directory
 *                 (as well as all the parent directorie created, if any).
 * @return An integer with 0 on success, last errno on failure
 */ 
int c_mkdir(const char *path, mode_t mode);


/**
 * Copy file to another.
 * @param to A C string with the new filepath
 * @param from A C string with the filepath to copy
 * @return An integer with 0 on success, 1 on failure.
 */
int c_copy(const char *to, const char *from);

/**
 * Remove file from filesytem
 * @note Also works for directory (if empty)
 * @param path A C string with the filepath to remove
 * @return An integer with 0 on success, last errno on failure
 */
int c_remove(const char *path);

/**
 * Remove a directory
 * @param path A C string with the path of the directory to remove.
 * @return An integer with 0 on success, last errno on failure
 */
int c_rmdir(const char *path);

/**
 * Remove a directory and its contents recursively
 *
 * This method mimics 'rm -rf' command. 
 * @param path A C string with the path of the directory to remove.
 * @return An integer with 0 on success, last errno on failure
 */
int c_rmdir_f(const char *path) ;

/**
 * Get some file informations
 * @note If the path cannot be "stat", most of the output parameters are set 
 * to -1.
 * @param[in] p A C string with the path of a file (or directory)
 * @param[out] pe An int with the permissions of the path 
 * @param[out] nl An int with the inumber of links
 * @param[out] ty An int with the type of the file :
 *    - 0 -> file
 *    - 1 -> link to a file
 *    - 2 -> directory
 *    - 3 -> link to a directory
 *    - 4 -> Other (fifo, socket, block special, char special ...)
 * @param[out] ui An int with the user id of the path 
 * @param[out] gi An int with the group id of the path 
 * @param[out] si An int with the size of the path 
 * @param[out] a A C string (20 chars wide, including NULL character) with the
 * last access date
 * @param[out] m A C string (20 chars wide, including NULL character) with the
 * last modification date
 * @param[out] c A C string (20 chars wide, including NULL character) with the 
 * creation date
 * @return An integer with 0 on success, last errno on failure 
 */
int c_fstat(const char *p, int *pe, int *nl, int *ty, int *ui, int *gi, 
            long *si, char a[20], char m[20], char c[20]);

/**
 * Check if path is accessible
 * @param[in] path A C string with the path to check
 * @param[in] perm An integer with the user's permission to check :
 *   - 0 do not check for permissions
 *   - 1 check for execute permission 
 *   - 2 check for write permission
 *   - 4 check for read permission
 * @return An int with 0 on success, errno on failure
 */
int c_access(const char *path, int perm) ;

/**
 * Create a directory or a file in given path
 * @param[in] path Path to be created
 * @param[in] mode Permission of the path
 * @param[in] astype A boolean with False to create a directory, True to create a file.
 * @param[in] forced A boolean with True to force creation of intermediate directory
 * @return 0 on success, -9 on allocation failure, last errno otherwise
 */
int c_create(const char* path, mode_t mode, int astype, int forced);

/**
 * Create a directory or a file in given path
 * @param[out] rows Number of rows of the current terminal window
 * @param[out] cols Number of columns of the current terminal window
 * @return An int with 0 on success, errno on failure. On failure, rows is set 
 * to 80 and cols to 20.
 */
int c_termsize(int *rows,int *cols);

/**
 * Get the current resident set size memory used by the program.
 */
size_t c_getCurrentRSS();

/**
 * Get the peak resident set size memory used by the program.
 */
size_t c_getPeakRSS();

/**
 * Get global memory usage informations.
 *
 * Note: The method attempts to read /proc/meminfo. If the file does not exists, all the given output arguments are
 * set to zero.
 */
int c_getSystemMemory(long long int *m_total,long long int *m_available,long long int *m_free);
