#!/usr/bin/awk -f
##########
#The MIT License (MIT)
#
# Copyright (c) 2015 Aiden Lab
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
##########

# Dedup script that submits deduping jobs after splitting at known 
# non-duplicate
# Juicer version 1.5
BEGIN{
    tot=0;
    name=0;
}
{
    if (tot >= 1000000) {
	if (p1 != $1 || p2 != $2 || p4 != $4 || p5 != $5 || p8 != $8) {
	    sname = sprintf("%s_msplit%04d_", groupname, name);
	    sscriptname = sprintf("%s/.%s.slurm", debugdir, sname);
	    if (justexact) {printf("awk -f %s/scripts/dups.awk -v nowobble=1 -v name=%s/%s %s/split%04d;\necho Reads:%s\n",  juicedir, dir, sname, dir, name, tot) > sscriptname;
        }
	    else {printf("awk -f %s/scripts/dups.awk -v name=%s/%s %s/split%04d; \necho Reads:%s \n", juicedir, dir, sname, dir, name, tot) > sscriptname;
	    }
	    sysstring = sprintf("cat %s", sscriptname);
	    system(sysstring);
	    outname = sprintf("%s/split%04d", dir, name);
	    close(outname);
	    close(sscriptname);
	    name++;
	    tot=0;
	}
    }
    outname = sprintf("%s/split%04d", dir, name);
    print > outname;
    p1=$1;p2=$2;p4=$4;p5=$5;p6=$6;p8=$8;
    tot++;
}
END {
    sname = sprintf("%s_msplit%04d_", groupname, name);
    sscriptname = sprintf("%s/.%s.slurm", debugdir, sname);
    if (justexact) {
	printf("awk -f %s/scripts/dups.awk -v nowobble=1 -v name=%s/%s %s/split%04d; \necho Reads:%s", juicedir, dir, sname, dir, name, tot) > sscriptname;
    }
    else {
	printf("awk -f %s/scripts/dups.awk -v name=%s/%s %s/split%04d; \necho Reads:%s \n", juicedir, dir, sname, dir, name, tot) > sscriptname;
    }
    sysstring = sprintf("cat %s", sscriptname);
    system(sysstring);
    close(sscriptname);
    sscriptname = sprintf("%s/.%s_msplit.slurm", debugdir, groupname);
    printf("\necho \"\"; cat %s/%s_msplit*_optdups.txt > %s/opt_dups.txt;cat %s/%s_msplit*_dups.txt > %s/dups.txt;cat %s/%s_msplit*_merged_nodups.txt > %s/merged_nodups.txt; \ndate \n", dir, groupname, dir, dir, groupname, dir, dir, groupname, dir) > sscriptname;
    sysstring = sprintf("cat %s", sscriptname);
    system(sysstring);
    close(sscriptname);
    sscriptname = sprintf("%s/.%s_finalize.slurm", debugdir, groupname);
    printf("echo %s %s %s %s;\n", topDir, site, genomeID, genomePath) > sscriptname;
    sysstring = sprintf("cat %s", sscriptname);
    system(sysstring);
    sysstring = sprintf("echo JobID=NONE dependency=afterany:NONE");
    system(sysstring);
    close(sscriptname);
    sscriptname = sprintf("%s/.%s_mail.slurm", debugdir, groupname);
    printf("echo %s %s %s %s\n", topDir, site, genomeID, genomePath) > sscriptname;
    sysstring = sprintf("cat %s", sscriptname);
    system(sysstring);
    close(sscriptname);
}

