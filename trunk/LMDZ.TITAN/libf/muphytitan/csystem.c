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

#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <libgen.h>
#include <sys/param.h>  // MAXPATHLEN
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/ioctl.h>
#include <pwd.h>
#include <grp.h>
#include <dirent.h>
#include <fcntl.h>
#include <time.h>
#include <unistd.h>
#include "csystem.h"

// Hack to get MAXPATHLEN
#ifndef MAXPATHLEN
#define MAXPATHLEN 4096
#endif

/* Get errno */
int c_get_errno(){
  int ret = errno;
  return ret;
}

/* Get the current working directory */
char* c_getcwd(){
  char *tmp;
  tmp = (char *) malloc((MAXPATHLEN+1)*sizeof(char));
  memset(tmp,'\0',MAXPATHLEN+1);
  if (getcwd(tmp,MAXPATHLEN) == NULL) {
    free(tmp);
    tmp = NULL;
  } 
  return tmp;
}


/* Get the realpath of input path and saves it in output path */ 
char* c_realpath(const char *input){
  if (!strlen(input)) {
    return c_getcwd();
  }
  return realpath(input,NULL);
}


/* Get directory name of input path */
char* c_dirname(const char *in){
  char *tmp = strdup(in);
  char *ou = strdup(dirname(tmp));
  free(tmp);
  return ou;
}

/* Get base name of input path */
char* c_basename(const char *in){
  char* tmp = strdup(in);
  char *ou = strdup(basename(tmp));
  free(tmp);
  return ou;
}

/* Get the corresponding name of the given user id */
char* c_uname(int uid){
  struct passwd  *pwd;
  char *name;
  name = NULL;
  if ((pwd = getpwuid(uid)) !=  NULL){
    name = strdup(pwd->pw_name);
  }
  return name;
}

/* Get the corresponding name of the given group id */
char* c_gname(int gid){
  struct group *grp;
  char *name;
  name = NULL;
  if ((grp = getgrgid(gid)) != NULL){
    name = strdup(grp->gr_name);
  }
  return name;
}

/* Get the error message of the given error id */
char* c_strerror(int err){
  char *tmp = strdup(strerror(err));
  if (errno != 0)
    tmp = "Unknown error\0";
  return tmp;
}


/* Remove a directory and its contents recursively */
int c_rmdir_f(const char *path) {
  DIR *d = opendir(path);
  size_t path_len = strlen(path);
  int r = -1;
  if (d) {
    struct dirent *p;
    r = 0;
    while (!r && (p=readdir(d))) {
      int r2 = -1;
      char *buf;
      size_t len;
      /* Skip the names "." and ".." as we don't want to recurse on them. */
      if (!strcmp(p->d_name, ".") || !strcmp(p->d_name, "..")) {
        continue;
      }
      len = path_len + strlen(p->d_name) + 2; 
      buf = malloc(len);
      if (buf) {
        struct stat statbuf;
        snprintf(buf, len, "%s/%s", path, p->d_name);
        if (!stat(buf, &statbuf)) {
          if (S_ISDIR(statbuf.st_mode)) {
            r2 = c_rmdir_f(buf);
          } else {
            r2 = unlink(buf);
          }
        }
        free(buf);
      }
      r = r2;
    }
    closedir(d);
  }
  if (!r) {
    r = rmdir(path);
	return r?errno:0;
  }
  return r;
}

/* Get some file informations */
int c_fstat(const char *p, int *pe, int *nl, int *ty, int *ui, int *gi, 
            long *si, char a[20], char m[20], char c[20]){
  struct stat stb;
  struct tm *t;
  int ret ;
  //char *tmp;
  char tmp[20];
  *pe = *ty = *ui = *gi = *si = -1 ;
  if (stat(p, &stb) != -1) {
    *pe = (int)(stb.st_mode & (S_IRWXU | S_IRWXG | S_IRWXO));
    *nl = (int)(stb.st_nlink);
    if (S_ISREG(stb.st_mode)) *ty = 0 + S_ISLNK(stb.st_mode);
    else if (S_ISDIR(stb.st_mode)) *ty = 2 + S_ISLNK(stb.st_mode) ;
    else *ty = 4;
    *ui = (int)stb.st_uid ; *gi = (int)stb.st_gid ; *si = (long)stb.st_size ;
    t = localtime(&stb.st_atime) ; ret = strftime(tmp, 20, "%F,%T", t); 
    if(ret != 0) {strncpy(a,tmp,20); a[19] = '\0';}else{a[0] = '\0';}
    t = localtime(&stb.st_mtime) ; ret = strftime(tmp, 20, "%F,%T", t);
    if(ret != 0) {strncpy(m,tmp,20); m[19] = '\0';}else{m[0] = '\0';}
    t = localtime(&stb.st_ctime) ; ret = strftime(tmp, 20, "%F,%T", t);
    if(ret != 0) {strncpy(c,tmp,20); c[19] = '\0';}else{c[0] = '\0';}
    return 0;
  }
  return errno;
}

/* Create a directory or a file in given path */
int c_create(const char* path, mode_t mode, int astype, int forced){
  char *p, *d;
  mode_t perm, cmask;
  int eval;
  // we build a directory
  if (!astype){
    if(forced){eval=c_mkdirp(path,mode);}else{eval=c_mkdir(path,mode);}
  }else{
    // we build a file
    if (forced){
      p = strdup(path) ; if(p == NULL) {return errno;}
      d = dirname(p) ; free(p) ; 
      if(d == NULL) {return -9;}
      // we attempts to create parent directory first
      cmask = umask(0);
      perm =((S_IRWXU | S_IRWXG | S_IRWXO) & ~(cmask & ~(S_IWUSR | S_IXUSR)));
      (void)umask(cmask) ;
      eval = c_mkdirp(d,perm); 
      if(eval){return eval;}
    }
    eval = open(path,O_CREAT|O_EXCL,mode);
    if (eval == -1) {eval=errno;}else{close(eval);eval=0;}
  }  
  return eval ;
}

/* Change current working directory */
int c_chdir(const char *path){
  if (chdir(path) != 0) return errno ;
  return 0;
}

/* Remove a directory */
int c_rmdir(const char *path){
 if ( rmdir(path) != 0) return errno;
  return 0;
}

/* Check if path is accessible */
int c_access(const char *path, int perm) {
  int ret;
  if (perm<0 && perm >7) perm=0;
  ret = access(path,perm);
  if (ret) return errno;
  return 0;
}

/* Rename a path */
int c_rename(const char *old, const char *new){
  if (rename(old,new) != 0) return errno ;
  return 0;
}

/* Change path permissions */
int c_chmod(const char *path, mode_t mode){
  if (chmod(path,mode) != 0) return errno ;
  return 0;
}

/* Create directory */
int c_mkdir(const char *path, mode_t mode){
  if(mkdir(path,mode) == 0) {
    if (chmod(path, mode) == -1) return errno ;
    return 0;
  }else{
    return errno;
  }
}


int c_copy(const char *to, const char *from) {
  int fd_to, fd_from;
  char buf[4096];
  ssize_t nread;
  int saved_errno;

  fd_from = open(from, O_RDONLY);
  if (fd_from < 0)
    return -1;

  fd_to = open(to, O_WRONLY | O_CREAT | O_EXCL, 0666);
  if (fd_to < 0)
    goto out_error;

  while (nread = read(fd_from, buf, sizeof buf), nread > 0) {
    char *out_ptr = buf;
    ssize_t nwritten;

    do {
      nwritten = write(fd_to, out_ptr, nread);

      if (nwritten >= 0) {
        nread -= nwritten;
        out_ptr += nwritten;
      } else if (errno != EINTR) {
        goto out_error;
      }
    } while (nread > 0);
  }

  if (nread == 0) {
    if (close(fd_to) < 0) {
      fd_to = -1;
      goto out_error;
    }
    close(fd_from);
    /* Success! */
    return 0;
  }

out_error:
  saved_errno = errno;

  close(fd_from);
  if (fd_to >= 0)
    close(fd_to);

  errno = saved_errno;
  return 1;
}

/* Remove file from filesytem */
int c_remove(const char *path){
  if (remove(path) != 0) return errno ;
  return 0;
}

/* Get terminal size */
int c_termsize(int *rows,int *cols){
  struct winsize max;
  int retval ;
  retval = ioctl(0, TIOCGWINSZ , &max);
  if (retval != 0) {
    retval = errno ;
    *rows = 20 ; *cols = 80 ;
  }else{
    *rows = max.ws_row ;
    *cols = max.ws_col ;
  }
  return retval;
}

/* Get the current umask (in decimal system) */
int c_umask(){
  mode_t mask ;
  mask = umask(S_IRWXU | S_IRWXG | S_IRWXO); (void) umask(mask) ;
  return (int)mask;
}

/* ------------------ THIRD PARTY ------------------ */

/*
 *  Copyright (c) 1983, 1992, 1993
 *   The Regents of the University of California.  All rights reserved.
 * 
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *  1. Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *  3. All advertising materials mentioning features or use of this software
 *     must display the following acknowledgement:
 *   This product includes software developed by the University of
 *  California, Berkeley and its contributors.
 *  4. Neither the name of the University nor the names of its contributors
 *     may be used to endorse or promote products derived from this software
 *     without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
 *  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 *  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 *  OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 *  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 *  OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 *  SUCH DAMAGE.
 */
/* Creates directories recursively */
int c_mkdirp(const char *path, mode_t mode) {
  struct stat sb;
  mode_t numask, oumask;
  int first, last, retval;
  char *cpath, *p;
  cpath=strdup(path); if (cpath == NULL) return -1;
  p = cpath; oumask = 0; retval = 0;
  if (p[0] == '/') ++p;    /* Skip leading '/'. */
  for (first = 1, last = 0; !last ; ++p) {
    if (p[0] == '\0') last = 1;
    else if (p[0] != '/') continue;
    *p = '\0';
    //if (! strlen(p) && p[1] == '\0') last = 1;
    if (first) {
      oumask = umask(0);
      numask = oumask & ~(S_IWUSR | S_IXUSR);
      (void)umask(numask);
      first = 0;
    }
    if (last) (void)umask(oumask);
    if (mkdir(cpath, last ? mode : S_IRWXU | S_IRWXG | S_IRWXO) < 0) {
      if (errno == EEXIST || errno == EISDIR) {
        if (stat(cpath, &sb) < 0) {
          retval = 1; break;
        } else if (!S_ISDIR(sb.st_mode)) {
          if (last) errno = EEXIST; else errno = ENOTDIR;
          retval = 1; break;
        }
      } else {
        retval = 1; break;
      }
    } 
    if (!last) *p = '/';
  }
  if (!first && !last) (void)umask(oumask);
  if (!retval){
    if (chmod(path, mode) == -1)  retval = errno ;
  } else {
   retval = errno;
  }
  free(cpath);
  return retval ;
}

/*
 * Copyright (c) 1999 Apple Computer, Inc. All rights reserved.
 *
 * @APPLE_LICENSE_HEADER_START@
 * 
 * "Portions Copyright (c) 1999 Apple Computer, Inc.  All Rights
 * Reserved.  This file contains Original Code and/or Modifications of
 * Original Code as defined in and that are subject to the Apple Public
 * Source License Version 1.0 (the 'License').  You may not use this file
 * except in compliance with the License.  Please obtain a copy of the
 * License at http://www.apple.com/publicsource and read it before using
 * this file.
 * 
 * The Original Code and all software distributed under the License are
 * distributed on an 'AS IS' basis, WITHOUT WARRANTY OF ANY KIND, EITHER
 * EXPRESS OR IMPLIED, AND APPLE HEREBY DISCLAIMS ALL SUCH WARRANTIES,
 * INCLUDING WITHOUT LIMITATION, ANY WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE OR NON-INFRINGEMENT.  Please see the
 * License for the specific language governing rights and limitations
 * under the License."
 * 
 * @APPLE_LICENSE_HEADER_END@
 *
 * Modified version: J. Burgalat (2015)
 *
 */

char* c_relpath(const char *path, const char *reldir) {
  char start_path[MAXPATHLEN+1];
  char end_path[MAXPATHLEN+1];
  char *buf;
  struct stat st;
  int i,mc;
  int last_elem;
  int prev_path;
  char *tmp,*cur;
  buf = NULL;

  if (realpath(reldir,start_path) == NULL) 
    return buf;
  
  if (realpath(path,end_path) == NULL) 
    return buf;

  // stat the starting path
  if (stat(start_path, &st) < 0) 
    return buf;
  
  if ((st.st_mode & S_IFMT) != S_IFDIR) {
    errno = ENOTDIR;
    return buf;
  } 
  if (start_path[strlen(start_path) - 1] != '/')
    strcat(start_path, "/");

  // stat the ending path path
  if (stat(end_path, &st) < 0) 
    return buf;
  
  if ((st.st_mode & S_IFMT) == S_IFDIR
      && end_path[strlen(end_path) - 1] != '/')
    strcat(end_path, "/");

  tmp = (char *) malloc((2*MAXPATHLEN+1)*sizeof(char));
  memset(tmp,'\0',2*MAXPATHLEN+1);
  cur = tmp;
  /* strip common prefix */
  i = 0;
  last_elem = 0;
  while (start_path[i] && start_path[i] == end_path[i]) {
    if (start_path[i] == '/')
      last_elem = i + 1;
    i += 1;
  }
  prev_path = 0;
  for (i = last_elem; start_path[i]; i += 1) {
    if (start_path[i] == '/') {
      if (prev_path){
        *cur = '/'; cur++;
      }
      strncpy(cur,"..",2) ; cur +=2;
      prev_path = 1;
    }
  }
  if (end_path[last_elem]) {
    if (prev_path) {
      *cur = '/'; cur++;
    }
    prev_path = 1;
    while (end_path[strlen(end_path) - 1] == '/')
      end_path[strlen(end_path) - 1] = '\0';
    strcpy(cur,&end_path[last_elem]) ;
    cur += strlen(&end_path[last_elem]);
  }
  if (! prev_path){
    *cur = '.' ; cur++ ;
  }
  // Normally here tmp contains our path. We just need to copy it in buf
  mc = strlen(tmp)+1;
  buf = (char *)malloc(mc*sizeof(char));
  memset(buf,'\0',mc);
  strncpy(buf,tmp,strlen(tmp));
  cur = buf;
  i=0;
  while (*cur != '\0' && i < mc){
    cur++;
    i++;
  }
  free(tmp);
  errno = 0;
  return buf;
}

/*
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 *
 * Cross-plateform way to get memory usage at a given time in fucking bytes !!
 *
 * Usage:
 *
 * size_t currentSize = getCurrentRSS( );
 * size_t peakSize    = getPeakRSS( );
 */

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <sys/resource.h>
#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>
#elif(defined(_AIX) || defined(__TOS__AIX__)) ||                                                                               \
    (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <procfs.h>
#endif
#endif

/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
size_t c_getPeakRSS() {
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
    return (size_t)info.PeakWorkingSetSize;

#elif(defined(_AIX) || defined(__TOS__AIX__)) ||                                                                               \
    (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
    /* AIX and Solaris ------------------------------------------ */
    struct psinfo psinfo;
    int fd = -1;
    if ((fd = open("/proc/self/psinfo", O_RDONLY)) == -1)
        return (size_t)0L; /* Can't open? */
    if (read(fd, &psinfo, sizeof(psinfo)) != sizeof(psinfo)) {
        close(fd);
        return (size_t)0L; /* Can't read? */
    }
    close(fd);
    return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    /* BSD, Linux, and OSX -------------------------------------- */
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
#if defined(__APPLE__) && defined(__MACH__)
    return (size_t)rusage.ru_maxrss;
#else
    return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
    /* Unknown OS ----------------------------------------------- */
    return (size_t)0L; /* Unsupported. */
#endif
}

/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
size_t c_getCurrentRSS() {
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
    return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
    /* OSX ------------------------------------------------------ */
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &infoCount) != KERN_SUCCESS)
        return (size_t)0L; /* Can't access? */
    return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    /* Linux ---------------------------------------------------- */
    long rss = 0L;
    FILE *fp = NULL;
    if ((fp = fopen("/proc/self/statm", "r")) == NULL)
        return (size_t)0L; /* Can't open? */
    if (fscanf(fp, "%*s%ld", &rss) != 1) {
        fclose(fp);
        return (size_t)0L; /* Can't read? */
    }
    fclose(fp);
    return (size_t)rss * (size_t)sysconf(_SC_PAGESIZE);

#else
    /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
    return (size_t)0L; /* Unsupported. */
#endif
}


 
/* Get some informatiosn about the OS memory usage (from /proc/meminfo) */
int c_getSystemMemory(long long int *m_total,long long int *m_available,long long int *m_free){
    FILE* fp;
    char buf[1024];
    char *tmp,*p;
    int ts;
    long long int mem1 = 0L,mem2 = 0L,mem3 = 0L;
    if (m_total) (*m_total) = mem1;
    if (m_available) (*m_free) = mem2;
    if (m_free) (*m_free) = mem3;
    if ((fp = fopen("/proc/meminfo", "r")) == NULL)  
        return 1;
    while (fgets(buf, sizeof(buf), fp) != NULL){
        ts = strlen(buf) - 1; buf[ts] = '\0'; 
        tmp = strndup(&buf[0],ts);
        // check our 3 cases
        if (ts >= 9 && !strncmp(tmp,"MemTotal:",9)){
            // extract the two first tokens.
            p=strtok(tmp, " "); p=strtok(NULL, " ");
            // convert the value.
            mem1 = strtoll(p,NULL,10);
        }
        if (ts >= 13 && !strncmp(tmp,"MemAvailable:",13)){
            p=strtok(tmp, " "); p=strtok(NULL, " "); 
            mem2 = strtoll(p,NULL,10);
        }    
        if (ts >= 8 && !strncmp(tmp,"MemFree:",8)){
            p=strtok(tmp, " "); p=strtok(NULL, " "); 
            mem3 = strtoll(p,NULL,10);
        }
        free(tmp);
    }
    fclose(fp);
    if (m_total) (*m_total) = mem1;
    if (m_available) (*m_available) = mem2;
    if (m_free) (*m_free) = mem3;
    return 0;
}
