from astropy.time import Time
from datetime import timedelta
from datetime import datetime
from datetime import date
import time
from astropy.table import Table
import numpy as np

import urllib
from astroplan import Observer
import astropy.units as u


class Weather:

    def daterange(self, start_date, end_date):
        for n in range(int((end_date - start_date).days)):
            yield start_date + timedelta(n)

    def read_davis_data_from_archive(self, date_obs, station="T60"):

        if station.capitalize == "T60":
            urlink = "http://t60meteo.tug.tubitak.gov.tr/index.html/"
            "Archive/ARC-%Y-%m-%d.txt".format
        elif station.capitalize == "RTT150":
            urlink = "http://rtt150meteo.tug.tubitak.gov.tr"
            "Archive/ARC-%Y-%m-%d.txt"
        elif station.capitalize == "T100":
            urlink = "http://t100meteo.tug.tubitak.gov.tr/index.html/"
            "Archive/ARC-%Y-%m-%d.txt"
        else:
            print("No station has found!")
            raise SystemExit

        raw_data = urllib.request.urlopen(urlink)
        dataset = np.genfromtxt(raw_data,
                                delimiter=None,
                                comments="---",
                                skip_header=1,
                                dtype="|U30")
        # date and time merging
        date_dataset = dataset.astype(object)
        date_dataset[:, 0] += "T" + date_dataset[:, 1]

        """
        convert first column time format because davis data
        not formatted properly
        """
        
        ts = [time.strftime('%Y-%m-%dT%H:%M',
                            time.strptime(s[0], '%Y%m%dT%H:%M')) for s in
              date_dataset[:, :1]]

        # convert np arrar for merging rest of data
        formatted_date_column = np.asarray(ts).reshape(len(ts), 1)

        np_davis_data = np.hstack((formatted_date_column,
                                   dataset[:, 2:]))

        tbl_davis_data = Table(np_davis_data,
                               names=('Timestamp',
                                      'Temp',
                                      'Chill',
                                      'HIindex',
                                      'Humid',
                                      'Dewpt',
                                      'Wind',
                                      'HiWind',
                                      'WindDir',
                                      'Rain',
                                      'Barom',
                                      'Solar',
                                      'ET',
                                      'UV'),
                               dtype=('datetime64',
                                      'f8',
                                      'f8',
                                      'f8',
                                      'f8',
                                      'f8',
                                      'f8',
                                      'f8',
                                      'f8',
                                      'f8',
                                      'f8',
                                      'f8',
                                      'f8',
                                      'f8'))
        
        return(tbl_davis_data)

    def astronomical_twilight(self, date_obs,
                              site_longitude=30.335555,
                              site_latitude=36.824166,
                              site_elevation=2500,
                              site_name="tug",
                              time_zone="Europe/Istanbul",
                              which="next"):

        # TUG's location info settings
        tug = Observer(longitude=site_longitude*u.deg,
                       latitude=site_latitude*u.deg,
                       elevation=site_elevation*u.m,
                       name=site_name,
                       timezone=time_zone)

        # convert date astropy date format
        astropy_time = Time(date_obs)

        # evening tw calculate
        et = tug.twilight_evening_astronomical(astropy_time,
                                               which=which)

        # morning tw calculate
        mt = tug.twilight_morning_astronomical(astropy_time,
                                               which=which)

        # localized time conversion
        return(tug.astropy_time_to_datetime(mt),
               tug.astropy_time_to_datetime(et))

    def calculate_bad_weather_time(self, date_obs, station="T60",
                                   site_longitude=30.335555,
                                   site_latitude=36.824166,
                                   site_elevation=2500,
                                   site_name="tug",
                                   time_zone="Europe/Istanbul"):

        # d = datetime.strptime(date_obs, "%Y-%m-%dT%H:%M:%S")
        # date_obs = d.strftime("%Y-%m-%d")

        mt, et = self.astronomical_twilight(date_obs,
                                            site_longitude=site_longitude,
                                            site_latitude=site_latitude,
                                            site_elevation=site_elevation,
                                            site_name=site_name,
                                            time_zone=time_zone)

        start_date = date(2017, 1, 1)
        end_date = date(2017, 12, 15)
