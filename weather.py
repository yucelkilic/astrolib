from datetime import timedelta
from datetime import date
import time

import numpy as np
import urllib
import re

from astropy.table import Table
from astropy.time import Time
import astropy.units as u
from astroplan import Observer


class Weather:

    def daterange(self, start_date, end_date):
        for n in range(int((end_date - start_date).days)):
            yield start_date + timedelta(n)

    def read_davis_data_from_archive(self, date_obs, station="T60"):
        """
        It reads wxview result data from meteo station.
        @param date_obs
        @type date_obs: date
        @param station: Station name: T60|T100|RTT150.
        @type station: string
        @return: astropy.table
        """

        year, month, day = re.split('[- :/.]', date_obs)

        if station.upper() == "T60":
            urlink = "http://t60meteo.tug.tubitak.gov.tr/index.html/" \
                     "Archive/ARC-{0}-{1}-{2}.txt".format(year, month, day)
        elif station.upper() == "RTT150":
            urlink = "http://rtt150meteo.tug.tubitak.gov.tr/" \
                     "ARC-{0}-{1}-{2}.txt".format(year, month, day)
        elif station.upper() == "T100":
            urlink = "http://t100meteo.tug.tubitak.gov.tr/index.html/" \
                     "Archive/ARC-{0}-{1}-{2}.txt".format(year, month, day)
        else:
            print("No station has found!")
            raise SystemExit

        # print("Retriving data from: {0}".format(urlink))
        try:
            raw_data = urllib.request.urlopen(urlink)
            dataset = np.genfromtxt(raw_data,
                                    delimiter=None,
                                    comments="---",
                                    skip_header=1,
                                    dtype="|U30")
        except:
            return False

        # date and time merging
        date_dataset = dataset.astype(object)
        date_dataset[:, 0] += "T" + date_dataset[:, 1]

        """
        convert first column time format because davis data
        not formatted properly
        """

        ts = [time.strftime('%Y-%m-%dT%H:%M',
                            time.strptime(s[0].replace("T24:00", "T23:59"),
                                          '%Y%m%dT%H:%M')) for s in
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

    def nautical_twilight(self, date_obs,
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
        et = tug.twilight_evening_nautical(astropy_time,
                                               which=which)

        # morning tw calculate
        mt = tug.twilight_morning_nautical(astropy_time,
                                               which=which)

        # localized time conversion
        return(tug.astropy_time_to_datetime(mt),
               tug.astropy_time_to_datetime(et))

    def daily_bad_weather_report(self, date_obs, station="T60",
                                 site_longitude=30.335555,
                                 site_latitude=36.824166,
                                 site_elevation=2500,
                                 site_name="tug",
                                 time_zone="Europe/Istanbul",
                                 humidity_limit=85,
                                 wind_limit=36):

        year, month, day = re.split('[- :/.]', date_obs)

        date_obs = "{0}-{1}-{2}".format(year, month, day)

        # d = datetime.strptime(date_obs, "%Y-%m-%dT%H:%M:%S")
        # date_obs = d.strftime("%Y-%m-%d")

        mt_before, et_before = self.nautical_twilight(
            date_obs,
            site_longitude,
            site_latitude,
            site_elevation,
            site_name,
            time_zone)

        davis_data_before_mid = self.read_davis_data_from_archive(
            date_obs,
            station=station)

        # date from after mid
        end_date = date(int(year), int(month), int(day)) + timedelta(1)
        mt_after, et_after = self.nautical_twilight(
            end_date.strftime("%Y-%m-%d"),
            site_longitude,
            site_latitude,
            site_elevation,
            site_name,
            time_zone)

        davis_data_after_mid = self.read_davis_data_from_archive(
            end_date.strftime("%Y-%m-%d"),
            station=station)

        # calculate dark night time
        tw_dates = [et_before, mt_after]
        t = Time(tw_dates)
        dark_hours = (t[1].jd - t[0].jd) * 24

        if davis_data_before_mid is not False:
            good_check = int((Time("{0}T23:59:00".format(date_obs)) - t[0]).sec / 300)
            coonection_check = int((Time("{0}T23:59:00".format(date_obs)) -
                                    Time(str(davis_data_before_mid['Timestamp'][len(davis_data_before_mid) - 1]))).sec / 300)
            if good_check > coonection_check:
                bad_data_before_mid = davis_data_before_mid[
                    ((davis_data_before_mid['Timestamp'] >= np.datetime64(et_before)) &
                     (davis_data_before_mid['Humid'] >= humidity_limit)) |
                    ((davis_data_before_mid['Timestamp'] >= np.datetime64(et_before)) &
                     (davis_data_before_mid['Wind'] >= wind_limit))]
                bad_data_before_mid = ["bad"] * (len(bad_data_before_mid) + coonection_check)
            else:
                bad_data_before_mid = ["bad"] * int((Time("{0}T23:59:00".format(date_obs)) - t[0]).sec / 300)
        else:
            bad_data_before_mid = ["bad"] *  int((Time("{0}T23:59:00".format(date_obs)) - t[0]).sec / 300)

        if davis_data_after_mid is not False:
            good_check = int((t[1] - Time("{0}T00:00:00".format(end_date.strftime("%Y-%m-%d")))).sec / 300)
            coonection_check = int((Time(str(davis_data_after_mid['Timestamp'][len(davis_data_after_mid) - 1])) -
                                    Time("{0}T00:00:00".format(end_date.strftime("%Y-%m-%d")))).sec / 300)
            if coonection_check > good_check:
                bad_data_after_mid = davis_data_after_mid[
                    ((davis_data_after_mid['Timestamp'] <= np.datetime64(mt_after)) &
                     (davis_data_after_mid['Humid'] >= humidity_limit)) |
                    ((davis_data_after_mid['Timestamp'] <= np.datetime64(mt_after)) &
                     (davis_data_after_mid['Wind'] >= wind_limit))]

                bad_data_after_mid = ["bad"] * (len(bad_data_after_mid))
            elif coonection_check < good_check:
                bad_data_after_mid = davis_data_after_mid[
                    ((davis_data_after_mid['Timestamp'] <= np.datetime64(mt_after)) &
                     (davis_data_after_mid['Humid'] >= humidity_limit)) |
                    ((davis_data_after_mid['Timestamp'] <= np.datetime64(mt_after)) &
                     (davis_data_after_mid['Wind'] >= wind_limit))]

                bad_data_after_mid = ["bad"] * (len(bad_data_after_mid) + coonection_check)
            else:
                bad_data_after_mid = ["bad"] * int(
                    (t[1] - Time("{0}T00:00:00".format(end_date.strftime("%Y-%m-%d")))).sec / 300)

        else:
            bad_data_after_mid = ["bad"] *  int(
                (t[1] - Time("{0}T00:00:00".format(end_date.strftime("%Y-%m-%d")))).sec / 300)

            # print(bad_data_before_mid['Timestamp', 'Humid', 'Wind'])
            # print(bad_data_after_mid['Timestamp', 'Humid', 'Wind'])

        total_bad_weather_time = ((len(bad_data_before_mid) +
                                   len(bad_data_after_mid)) * 5) / 60.0

        return(date_obs,
               '{0:.2f}'.format(dark_hours),
               '{0:.2f}'.format(total_bad_weather_time))

    def long_term_bad_weather_report(self, start_date,
                                     end_date,
                                     station="T60",
                                     site_longitude=30.335555,
                                     site_latitude=36.824166,
                                     site_elevation=2500,
                                     site_name="tug",
                                     time_zone="Europe/Istanbul",
                                     humidity_limit=85,
                                     wind_limit=36):
        syear, smonth, sday = re.split('[- :/.]', start_date)

        start_date = "{0}-{1}-{2}".format(syear, smonth, sday)

        eyear, emonth, eday = re.split('[- :/.]', end_date)

        end_date = "{0}-{1}-{2}".format(eyear, emonth, eday)

        start_date = date(int(syear), int(smonth), int(sday))
        end_date = date(int(eyear), int(emonth), int(eday))

        print("{0} {1} {2}\n".format("date", "dark_hours",
                                     "bad_weather_hours"))

        with open("{0}_bad_weather_report-{1}_{2}.txt".format(
                station,
                start_date,
                end_date), "a") as out:

            out.write("{0} {1} {2}\n".format("date", "dark_hours",
                                             "bad_weather_hours"))

            for single_date in self.daterange(start_date, end_date):
                odate = single_date.strftime("%Y-%m-%d")
                dt, dark_hours, bad_time = self.daily_bad_weather_report(
                    odate,
                    station,
                    site_longitude,
                    site_latitude,
                    site_elevation,
                    site_name,
                    time_zone,
                    humidity_limit,
                    wind_limit)

                print(dt, dark_hours, bad_time)
                out.write("{0} {1} {2}\n".format(dt,
                                                 dark_hours,
                                                 bad_time))

        print("Long term bad weather report has been written "
              "to {0}_bad_weather_report-{1}_{2}.txt".format(
                  station,
                  start_date,
                  end_date))

        return(True)
