namespace gevolution{


// writes lightcones of particles
void writeLightcones (
    metadata &sim, const cosmology cosmo, const double a,
    const double tau, const double dtau, const double dtau_old,
    const double /* maxvel */, 
    const int cycle, 
    std::string h5filename,
    Particles_gevolution *pcls_cdm,
    
    // redundant particle structures
    // Particles_gevolution *pcls_b,
    // Particles_gevolution *pcls_ncdm,
    
    // lightcone function for particles only
    // Field<Real> *phi, 
    // Field<Real> * /* chi */, 
    // Field<Real> * /* Bi */, 
    // Field<Real> * /* Sij */,
    // Field<Cplx> * /* BiFT */, 
    // Field<Cplx> * /* SijFT */, 
    // PlanFFT<Cplx> * /* plan_Bi */,
    // PlanFFT<Cplx> * /* plan_Sij */, 
    
    int &done_hij, 
    std::set<long> *IDbacklog)
{
    int i, j, n, p;
    double d;
    double vertex[MAX_INTERSECTS][3];
    double domain[6];
    // double pos[3];
    double s[2];
    char filename[2 * PARAM_MAX_LENGTH + 24];
    char buffer[268];
    FILE *outfile = nullptr;
    gadget2_header hdr;
    std::set<long> IDprelog[MAX_PCL_SPECIES];
    long *IDcombuf;
    long *IDcombuf2;
    Site xsim;
    // int done_B = 0;


    done_hij = 0;

    domain[0] = -0.5;
    domain[1] = phi->lattice ().coordSkip ()[1] - 0.5;
    domain[2] = phi->lattice ().coordSkip ()[0] - 0.5;
    for (i = 0; i < 3; i++)
        domain[i + 3] = domain[i] + phi->lattice ().sizeLocal (i) + 1.;

    for (i = 0; i < 6; i++)
        domain[i] /= (double)sim.numpts;

    for (i = 0; i < sim.num_lightcone; i++)
    {
        if (parallel.isRoot ())
        {
            if (sim.num_lightcone > 1)
                sprintf (filename, "%s%s%d_info.dat", sim.output_path,
                         sim.basename_lightcone, i);
            else
                sprintf (filename, "%s%s_info.dat", sim.output_path,
                         sim.basename_lightcone);

            outfile = fopen (filename, "a");
            if (outfile == NULL)
            {
                std::cout << " error opening file for lightcone info!"
                          << std::endl;
            }
            else if (cycle == 0)
            {
                if (sim.num_lightcone > 1)
                    fprintf (outfile,
                             "# information file for lightcone %d\n# "
                             "geometric "
                             "parameters:\n# vertex = (%f, %f, %f) "
                             "Mpc/h\n# "
                             "redshift = "
                             "%f\n# distance = (%f - %f) Mpc/h\n# "
                             "opening "
                             "half-angle = "
                             "%f degrees\n# direction = (%f, %f, %f)\n# "
                             "cycle   "
                             "tau/boxsize    a              pcl_inner   "
                             "     "
                             "pcl_outer  "
                             "      metric_inner     metric_outer\n",
                             i, sim.lightcone[i].vertex[0] * sim.boxsize,
                             sim.lightcone[i].vertex[1] * sim.boxsize,
                             sim.lightcone[i].vertex[2] * sim.boxsize,
                             sim.lightcone[i].z,
                             sim.lightcone[i].distance[0] * sim.boxsize,
                             sim.lightcone[i].distance[1] * sim.boxsize,
                             (sim.lightcone[i].opening > -1.)
                                 ? acos (sim.lightcone[i].opening) * 180. / M_PI
                                 : 180.,
                             sim.lightcone[i].direction[0],
                             sim.lightcone[i].direction[1],
                             sim.lightcone[i].direction[2]);
                else
                    fprintf (outfile,
                             "# information file for lightcone\n# "
                             "geometric "
                             "parameters:\n# vertex = (%f, %f, %f) "
                             "Mpc/h\n# "
                             "redshift = "
                             "%f\n# distance = (%f - %f) Mpc/h\n# "
                             "opening "
                             "half-angle = "
                             "%f degrees\n# direction = (%f, %f, %f)\n# "
                             "cycle   "
                             "tau/boxsize    a              pcl_inner   "
                             "     "
                             "pcl_outer  "
                             "      metric_inner     metric_outer\n",
                             sim.lightcone[i].vertex[0] * sim.boxsize,
                             sim.lightcone[i].vertex[1] * sim.boxsize,
                             sim.lightcone[i].vertex[2] * sim.boxsize,
                             sim.lightcone[i].z,
                             sim.lightcone[i].distance[0] * sim.boxsize,
                             sim.lightcone[i].distance[1] * sim.boxsize,
                             (sim.lightcone[i].opening > -1.)
                                 ? acos (sim.lightcone[i].opening) * 180. / M_PI
                                 : 180.,
                             sim.lightcone[i].direction[0],
                             sim.lightcone[i].direction[1],
                             sim.lightcone[i].direction[2]);
            }
        }

        d = particleHorizon (1. / (1. + sim.lightcone[i].z), cosmo);

        s[0] = d - tau - 0.5 * sim.covering[i] * dtau;
        s[1] = d - tau + 0.5 * sim.covering[i] * dtau_old;


        if (sim.lightcone[i].distance[0] > s[0]
            && sim.lightcone[i].distance[1] <= s[1] && s[1] > 0.)
        {
            if (parallel.isRoot () && outfile != NULL)
            {
                fprintf (outfile,
                         "%6d   %e   %e   %2.12f   %2.12f   %2.12f "
                         "  %2.12f\n",
                         cycle, tau, a, d - tau - 0.5 * dtau,
                         d - tau + 0.5 * dtau_old, s[0], s[1]);
                fclose (outfile);

                if (sim.num_lightcone > 1)
                    sprintf (filename, "%s%s%d_info.bin", sim.output_path,
                             sim.basename_lightcone, i);
                else
                    sprintf (filename, "%s%s_info.bin", sim.output_path,
                             sim.basename_lightcone);

                outfile = fopen (filename, "a");
                if (outfile == NULL)
                {
                    std::cout << " error opening file for lightcone "
                                 "info!"
                              << std::endl;
                }
                else
                {
                    ((double *)buffer)[0] = tau;
                    ((double *)buffer)[1] = a;
                    ((double *)buffer)[2] = d - tau - 0.5 * dtau;
                    ((double *)buffer)[3] = d - tau + 0.5 * dtau_old;

                    fwrite ((const void *)&cycle, sizeof (int), 1, outfile);
                    fwrite ((const void *)buffer, sizeof (double), 4, outfile);
                    fwrite ((const void *)s, sizeof (double), 2, outfile);

                    fclose (outfile);
                }
            }

        }
        else if (parallel.isRoot () && outfile != NULL)
            fclose (outfile);

        if (sim.out_lightcone[i] & MASK_GADGET
            && sim.lightcone[i].distance[0] > d - tau + 0.5 * dtau_old
            && sim.lightcone[i].distance[1] <= d - tau + 0.5 * dtau_old
            && d - tau + 0.5 * dtau_old > 0.)
        {
            n = findIntersectingLightcones (
                sim.lightcone[i],
                d - tau + (0.5 + LIGHTCONE_IDCHECK_ZONE) * dtau_old,
                d - tau - 0.5 * dtau, domain, vertex);

            hdr.num_files = 1;
            hdr.Omega0 = cosmo.Omega_m;
            hdr.OmegaLambda = cosmo.Omega_Lambda;
            hdr.HubbleParam = cosmo.h;
            hdr.BoxSize = sim.boxsize / GADGET_LENGTH_CONVERSION;
            hdr.flag_sfr = 0;
            hdr.flag_cooling = 0;
            hdr.flag_feedback = 0;
            hdr.flag_age = 0;
            hdr.flag_metals = 0;
            for (p = 0; p < 256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4
                                - 4 * 8 - 2 * 4 - 6 * 4;
                 p++)
                hdr.fill[p] = 0;
            for (p = 0; p < 6; p++)
            {
                hdr.npart[p] = 0;
                hdr.npartTotal[p] = 0;
                hdr.npartTotalHW[p] = 0;
                hdr.mass[p] = 0.;
            }

            hdr.time = a;
            hdr.redshift = (1. / a) - 1.;

            if (sim.baryon_flag)
                hdr.mass[1] = (double)sim.tracer_factor[0] * cosmo.C_RHO_CRIT
                              * cosmo.Omega_cdm * sim.boxsize * sim.boxsize
                              * sim.boxsize / sim.numpcl[0]
                              / GADGET_MASS_CONVERSION;
            else
                hdr.mass[1] = (double)sim.tracer_factor[0] * cosmo.C_RHO_CRIT
                              * (cosmo.Omega_cdm + cosmo.Omega_b) * sim.boxsize
                              * sim.boxsize * sim.boxsize / sim.numpcl[0]
                              / GADGET_MASS_CONVERSION;

            if (sim.num_lightcone > 1)
                sprintf (filename, "%d_%04d", i, cycle);
            else
                sprintf (filename, "_%04d", cycle);

            if (sim.tracer_factor[0] > 0)
                pcls_cdm->saveGadget2 (h5filename + filename + "_cdm", hdr,
                                       sim.lightcone[i], d - tau, dtau,
                                       dtau_old, a * Hconf (a, cosmo),
                                       vertex, n, IDbacklog[0], IDprelog[0],
                                       phi, sim.tracer_factor[0]);

            if (sim.baryon_flag && sim.tracer_factor[1] > 0)
            {
                hdr.mass[1] = (double)sim.tracer_factor[1] * cosmo.C_RHO_CRIT
                              * cosmo.Omega_b * sim.boxsize * sim.boxsize
                              * sim.boxsize / sim.numpcl[1]
                              / GADGET_MASS_CONVERSION;
                pcls_b->saveGadget2 (h5filename + filename + "_b", hdr,
                                     sim.lightcone[i], d - tau, dtau, dtau_old,
                                     a * Hconf (a, cosmo), vertex, n,
                                     IDbacklog[1], IDprelog[1], phi,
                                     sim.tracer_factor[1]);
            }

            for (p = 0; p < cosmo.num_ncdm; p++)
            {
                if (sim.numpcl[1 + sim.baryon_flag + p] == 0
                    || sim.tracer_factor[p + 1 + sim.baryon_flag] == 0)
                    continue;
                sprintf (buffer, "_ncdm%d", p);
                hdr.mass[1] = (double)sim.tracer_factor[p + 1 + sim.baryon_flag]
                              * cosmo.C_RHO_CRIT * cosmo.Omega_ncdm[p] * sim.boxsize
                              * sim.boxsize * sim.boxsize
                              / sim.numpcl[p + 1 + sim.baryon_flag]
                              / GADGET_MASS_CONVERSION;
                pcls_ncdm[p].saveGadget2 (
                    h5filename + filename + buffer, hdr, sim.lightcone[i],
                    d - tau, dtau, dtau_old, a * Hconf (a, cosmo),
                    vertex, n, IDbacklog[p + 1 + sim.baryon_flag],
                    IDprelog[p + 1 + sim.baryon_flag], phi,
                    sim.tracer_factor[p + 1 + sim.baryon_flag]);
            }
        }
    }


    for (p = 0; p <= cosmo.num_ncdm + sim.baryon_flag; p++)
    {
        IDbacklog[p] = IDprelog[p];
        IDprelog[p].clear ();

        n = IDbacklog[p].size ();
        // dim 0 send/rec
        if (parallel.grid_rank ()[0] % 2 == 0)
        {
            parallel.send_dim0<int> (
                n, (parallel.grid_size ()[0] + parallel.grid_rank ()[0] - 1)
                       % parallel.grid_size ()[0]);
            parallel.receive_dim0<int> (
                i, (parallel.grid_size ()[0] + parallel.grid_rank ()[0] - 1)
                       % parallel.grid_size ()[0]);
            parallel.send_dim0<int> (n, (parallel.grid_rank ()[0] + 1)
                                            % parallel.grid_size ()[0]);
            parallel.receive_dim0<int> (j, (parallel.grid_rank ()[0] + 1)
                                               % parallel.grid_size ()[0]);
        }
        else
        {
            parallel.receive_dim0<int> (i, (parallel.grid_rank ()[0] + 1)
                                               % parallel.grid_size ()[0]);
            parallel.send_dim0<int> (n, (parallel.grid_rank ()[0] + 1)
                                            % parallel.grid_size ()[0]);
            parallel.receive_dim0<int> (
                j, (parallel.grid_size ()[0] + parallel.grid_rank ()[0] - 1)
                       % parallel.grid_size ()[0]);
            parallel.send_dim0<int> (
                n, (parallel.grid_size ()[0] + parallel.grid_rank ()[0] - 1)
                       % parallel.grid_size ()[0]);
        }

        if (n + i + j > 0)
        {
            IDcombuf = (long *)malloc ((n + i + j) * sizeof (long));

            n = 0;
            for (std::set<long>::iterator it = IDbacklog[p].begin ();
                 it != IDbacklog[p].end (); it++)
                IDcombuf[n++] = *it;

            if (parallel.grid_rank ()[0] % 2 == 0)
            {
                if (n > 0)
                    parallel.send_dim0<long> (IDcombuf, n,
                                              (parallel.grid_size ()[0]
                                               + parallel.grid_rank ()[0] - 1)
                                                  % parallel.grid_size ()[0]);
                if (i > 0)
                    parallel.receive_dim0<long> (
                        IDcombuf + n, i,
                        (parallel.grid_size ()[0] + parallel.grid_rank ()[0]
                         - 1)
                            % parallel.grid_size ()[0]);
                if (n > 0)
                    parallel.send_dim0<long> (IDcombuf, n,
                                              (parallel.grid_rank ()[0] + 1)
                                                  % parallel.grid_size ()[0]);
                if (j > 0)
                    parallel.receive_dim0<long> (
                        IDcombuf + n + i, j,
                        (parallel.grid_rank ()[0] + 1)
                            % parallel.grid_size ()[0]);
            }
            else
            {
                if (i > 0)
                    parallel.receive_dim0<long> (
                        IDcombuf + n, i,
                        (parallel.grid_rank ()[0] + 1)
                            % parallel.grid_size ()[0]);
                if (n > 0)
                    parallel.send_dim0<long> (IDcombuf, n,
                                              (parallel.grid_rank ()[0] + 1)
                                                  % parallel.grid_size ()[0]);
                if (j > 0)
                    parallel.receive_dim0<long> (
                        IDcombuf + n + i, j,
                        (parallel.grid_size ()[0] + parallel.grid_rank ()[0]
                         - 1)
                            % parallel.grid_size ()[0]);
                if (n > 0)
                    parallel.send_dim0<long> (IDcombuf, n,
                                              (parallel.grid_size ()[0]
                                               + parallel.grid_rank ()[0] - 1)
                                                  % parallel.grid_size ()[0]);
            }

            n += i + j;

            for (i = IDbacklog[p].size (); i < n; i++)
                IDbacklog[p].insert (IDcombuf[i]);
        }

        // dim 1 send/rec
        if (parallel.grid_rank ()[1] % 2 == 0)
        {
            parallel.send_dim1<int> (
                n, (parallel.grid_size ()[1] + parallel.grid_rank ()[1] - 1)
                       % parallel.grid_size ()[1]);
            parallel.receive_dim1<int> (
                i, (parallel.grid_size ()[1] + parallel.grid_rank ()[1] - 1)
                       % parallel.grid_size ()[1]);
            parallel.send_dim1<int> (n, (parallel.grid_rank ()[1] + 1)
                                            % parallel.grid_size ()[1]);
            parallel.receive_dim1<int> (j, (parallel.grid_rank ()[1] + 1)
                                               % parallel.grid_size ()[1]);

            if (n > 0)
                parallel.send_dim1<long> (
                    IDcombuf, n,
                    (parallel.grid_size ()[1] + parallel.grid_rank ()[1] - 1)
                        % parallel.grid_size ()[1]);

            if (i > 0)
            {
                IDcombuf2 = (long *)malloc (i * sizeof (long));
                parallel.receive_dim1<long> (
                    IDcombuf2, i,
                    (parallel.grid_size ()[1] + parallel.grid_rank ()[1] - 1)
                        % parallel.grid_size ()[1]);
                while (i > 0)
                    IDbacklog[p].insert (IDcombuf2[--i]);
                free (IDcombuf2);
            }

            if (n > 0)
            {
                parallel.send_dim1<long> (IDcombuf, n,
                                          (parallel.grid_rank ()[1] + 1)
                                              % parallel.grid_size ()[1]);
                free (IDcombuf);
            }

            if (j > 0)
            {
                IDcombuf2 = (long *)malloc (j * sizeof (long));
                parallel.receive_dim1<long> (IDcombuf2, j,
                                             (parallel.grid_rank ()[1] + 1)
                                                 % parallel.grid_size ()[1]);
                while (j > 0)
                    IDbacklog[p].insert (IDcombuf2[--j]);
                free (IDcombuf2);
            }
        }
        else
        {
            parallel.receive_dim1<int> (i, (parallel.grid_rank ()[1] + 1)
                                               % parallel.grid_size ()[1]);
            parallel.send_dim1<int> (n, (parallel.grid_rank ()[1] + 1)
                                            % parallel.grid_size ()[1]);
            parallel.receive_dim1<int> (
                j, (parallel.grid_size ()[1] + parallel.grid_rank ()[1] - 1)
                       % parallel.grid_size ()[1]);
            parallel.send_dim1<int> (
                n, (parallel.grid_size ()[1] + parallel.grid_rank ()[1] - 1)
                       % parallel.grid_size ()[1]);

            if (i > 0)
            {
                IDcombuf2 = (long *)malloc (i * sizeof (long));
                parallel.receive_dim1<long> (IDcombuf2, i,
                                             (parallel.grid_rank ()[1] + 1)
                                                 % parallel.grid_size ()[1]);
                while (i > 0)
                    IDbacklog[p].insert (IDcombuf2[--i]);
                free (IDcombuf2);
            }

            if (n > 0)
                parallel.send_dim1<long> (IDcombuf, n,
                                          (parallel.grid_rank ()[1] + 1)
                                              % parallel.grid_size ()[1]);

            if (j > 0)
            {
                IDcombuf2 = (long *)malloc (j * sizeof (long));
                parallel.receive_dim1<long> (
                    IDcombuf2, j,
                    (parallel.grid_size ()[1] + parallel.grid_rank ()[1] - 1)
                        % parallel.grid_size ()[1]);
                while (j > 0)
                    IDbacklog[p].insert (IDcombuf2[--j]);
                free (IDcombuf2);
            }

            if (n > 0)
            {
                parallel.send_dim1<long> (
                    IDcombuf, n,
                    (parallel.grid_size ()[1] + parallel.grid_rank ()[1] - 1)
                        % parallel.grid_size ()[1]);
                free (IDcombuf);
            }
        }
    }
}


} // namespace gevolution
