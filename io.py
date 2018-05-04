# -*- coding: utf-8 -*-

import glob
import numpy as np
import paramiko
import os
import sqlite3
from .astronomy import FitsOps
from datetime import datetime


class FileOps:

    def make_date(self, datestr, datefrmt='%Y-%m-%d'):

        return(datetime.strptime(datestr.decode('ascii'), datefrmt))

    def read_file_as_array(self, file_name):

        """
        Reads text file into numpy array.
        @param file_name: Text file name and path
        @type file_name: str
        @return: array
        """

        try:
            return (np.genfromtxt(file_name,
                                  comments='#',
                                  delimiter=' | ',
                                  dtype="U"))
        except Exception as e:
            print(e)

    def read_res(self, file_name):

        """
        Reads A-Track result file into numpy array.
        @param file_name: Text file name and path
        @type file_name: str
        @return: array
        """

        try:
            data = np.genfromtxt(file_name,
                                 comments='#',
                                 skip_header=2,
                                 invalid_raise=False,
                                 delimiter=None,
                                 usecols=(0, 1, 3, 4, 5))
            return (data[~np.isnan(data).any(axis=1)])
        except Exception as e:
            print(e)

    def get_file_list(self, dir_name):

        """
        List FITS images in a folder into a numpy array.
        @param dir_name: Directory of FITS images
        @type dir_name: str
        @return: array
        """

        images = sorted(glob.glob(dir_name + '/*.fit*'))
        return (images)

    def get_fits_from_server(self,
                             hostname,
                             username,
                             password,
                             dirname="/mnt/data/images",
                             fits_ext=".fts",
                             header2sqlite=False,
                             sqlite_file="gozlemler.sqlite"):
        
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

        try:
            ssh.connect(hostname=hostname,
                        username=username,
                        password=password)
        except paramiko.SSHException:
            print("Connection Failed")
            quit()

        sftp = ssh.open_sftp()
        sftp.chdir(dirname)

        ret = False

        for fileattr in sftp.listdir_attr():
            if not os.path.exists(fileattr.filename) and \
               fits_ext in fileattr.filename:
                sftp.get(fileattr.filename, fileattr.filename)
                if header2sqlite is True:
                    self.fitshead_to_database(fileattr.filename, sqlite_file=sqlite_file)
                ret = True
                print("{0} => {1}".format(fileattr.filename,
                                          fileattr.filename))

        print("Done")
        ssh.close()
        if ret is False:
            print("No file(s) found!")

        return(ret)

    def fitshead_to_database(self, fits_file,
                             sqlite_file="gozlemler.sqlite",
                             table_name="gozlemler",
                             keywords=None):

        if keywords is None:
            keywords = ['xfactor',
                        'yfactor',
                        'exptime',
                        'object',
                        'priority',
                        'instrume',
                        'jd',
                        'date-obs',
                        'time-obs',
                        'lst',
                        'latitude',
                        'elevatio',
                        'azimuth',
                        'ha',
                        'ra',
                        'dec',
                        'objra',
                        'objdec',
                        'epoch',
                        'equinox',
                        'filter',
                        'camtemp',
                        'focuspos',
                        'wxtemp',
                        'wxpres',
                        'wxwndspd',
                        'wxwnddir',
                        'wxhumid',
                        'biascor',
                        'thermcor',
                        'flatcor',
                        'badpxcor',
                        'fwhmh',
                        'fwhmhs',
                        'fwhmv',
                        'fwhmvs']

        # Connecting to the database file
        conn = sqlite3.connect(sqlite_file)
        c = conn.cursor()

        fo = FitsOps(fits_file)
        keyword_values = []
        for keyword in keywords:
            keyword_value = fo.get_header(keyword)
            if keyword_value is None:
                keyword_value = -9999
                    
            keyword_values.append(keyword_value)
                            
        c.execute("INSERT OR IGNORE INTO {tn} {cn} VALUES {vls}".format(
            tn=table_name,
            cn=tuple(keywords),
            vls=tuple(keyword_values)))

        conn.commit()
        conn.close()

        return(True)

    def find_if_in_database_id(self, database, idd):
        """
        Search detected asteroids ID in the MPCORB.DAT database for MPC report.
        @param database: MPCORB.DAT path
        @type database: str
        @param idd: Asteroid's ID
        @type idd: str
        @return: str
        """

        ret = ""
        try:
            f = open(database, "r")
            for i in f:
                ln = i.replace("\n", "").split()
                try:
                    if "({0})".format(idd) == ln[21]:
                        id_name = ln[0]
                        if len(id_name) > 5:
                            ret = "     " + id_name
                        else:
                            ret = id_name
                except:
                    continue
            f.close()
        except Exception as e:
            print(e)

        return (ret)

    def find_if_in_database_name(self, database, name):

        """
        Search detected asteroids ID by the name in the
        MPCORB.DAT database for MPC report.

        @param database: MPCORB.DAT path
        @type database: str
        @param name: Asteroid's name
        @type name: str
        @return: str
        """

        ret = ""
        try:
            f = open(database, "r")
            for i in f:
                ln = i.replace("\n", "").split()
                try:
                    if len(ln[23]) < 8:
                        combname = "{0} {1}".format(ln[22], ln[23])
                    else:
                        combname = ln[22]

                    if name == combname:
                        id_name = ln[0]
                        if len(id_name) > 5:
                            ret = "     " + id_name
                        else:
                            ret = id_name
                except:
                    continue
            f.close()
        except Exception as e:
            print(e)

        return (ret)
