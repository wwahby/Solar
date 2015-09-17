function [flux_per_m2] = get_photon_flux_per_m2( lambda_nm )
    lambda_micron_vec = [0.1195	0.1205	0.1215	0.1225	0.1235	0.1245	0.1255	0.1265	0.1275	0.1285	0.1295	0.1305	0.1315	0.1325	0.1335	0.1345	0.1355	0.1365	0.1375	0.1385	0.1395	0.1405	0.1415	0.1425	0.1435	0.1445	0.1455	0.1465	0.1475	0.1485	0.1495	0.1505	0.1515	0.1525	0.1535	0.1545	0.1555	0.1565	0.1575	0.1585	0.1595	0.1605	0.1615	0.1625	0.1635	0.1645	0.1655	0.1665	0.1675	0.1685	0.1695	0.1705	0.1715	0.1725	0.1735	0.1745	0.1755	0.1765	0.1775	0.1785	0.1795	0.1805	0.1815	0.1825	0.1835	0.1845	0.1855	0.1865	0.1875	0.1885	0.1895	0.1905	0.1915	0.1925	0.1935	0.1945	0.1955	0.1965	0.1975	0.1985	0.1995	0.2005	0.2015	0.2025	0.2035	0.2045	0.2055	0.2065	0.2075	0.2085	0.2095	0.2105	0.2115	0.2125	0.2135	0.2145	0.2155	0.2165	0.2175	0.2185	0.2195	0.2205	0.2215	0.2225	0.2235	0.2245	0.2255	0.2265	0.2275	0.2285	0.2295	0.2305	0.2315	0.2325	0.2335	0.2345	0.2355	0.2365	0.2375	0.2385	0.2395	0.2405	0.2415	0.2425	0.2435	0.2445	0.2455	0.2465	0.2475	0.2485	0.2495	0.2505	0.2515	0.2525	0.2535	0.2545	0.2555	0.2565	0.2575	0.2585	0.2595	0.2605	0.2615	0.2625	0.2635	0.2645	0.2655	0.2665	0.2675	0.2685	0.2695	0.2705	0.2715	0.2725	0.2735	0.2745	0.2755	0.2765	0.2775	0.2785	0.2795	0.2805	0.2815	0.2825	0.2835	0.2845	0.2855	0.2865	0.2875	0.2885	0.2895	0.2905	0.2915	0.2925	0.2935	0.2945	0.2955	0.2965	0.2975	0.2985	0.2995	0.3005	0.3015	0.3025	0.3035	0.3045	0.3055	0.3065	0.3075	0.3085	0.3095	0.3105	0.3115	0.3125	0.3135	0.3145	0.3155	0.3165	0.3175	0.3185	0.3195	0.3205	0.3215	0.3225	0.3235	0.3245	0.3255	0.3265	0.3275	0.3285	0.3295	0.3305	0.3315	0.3325	0.3335	0.3345	0.3355	0.3365	0.3375	0.3385	0.3395	0.3405	0.3415	0.3425	0.3435	0.3445	0.3455	0.3465	0.3475	0.3485	0.3495	0.3505	0.3515	0.3525	0.3535	0.3545	0.3555	0.3565	0.3575	0.3585	0.3595	0.3605	0.3615	0.3625	0.3635	0.3645	0.3655	0.3665	0.3675	0.3685	0.3695	0.3705	0.3715	0.3725	0.3735	0.3745	0.3755	0.3765	0.3775	0.3785	0.3795	0.3805	0.3815	0.3825	0.3835	0.3845	0.3855	0.3865	0.3875	0.3885	0.3895	0.3905	0.3915	0.3925	0.3935	0.3945	0.3955	0.3965	0.3975	0.3985	0.3995	0.4005	0.4015	0.4025	0.4035	0.4045	0.4055	0.4065	0.4075	0.4085	0.4095	0.4105	0.4115	0.4125	0.4135	0.4145	0.4155	0.4165	0.4175	0.4185	0.4195	0.4205	0.4215	0.4225	0.4235	0.4245	0.4255	0.4265	0.4275	0.4285	0.4295	0.4305	0.4315	0.4325	0.4335	0.4345	0.4355	0.4365	0.4375	0.4385	0.4395	0.4405	0.4415	0.4425	0.4435	0.4445	0.4455	0.4465	0.4475	0.4485	0.4495	0.4505	0.4515	0.4525	0.4535	0.4545	0.4555	0.4565	0.4575	0.4585	0.4595	0.4605	0.4615	0.4625	0.4635	0.4645	0.4655	0.4665	0.4675	0.4685	0.4695	0.4705	0.4715	0.4725	0.4735	0.4745	0.4755	0.4765	0.4775	0.4785	0.4795	0.4805	0.4815	0.4825	0.4835	0.4845	0.4855	0.4865	0.4875	0.4885	0.4895	0.4905	0.4915	0.4925	0.4935	0.4945	0.4955	0.4965	0.4975	0.4985	0.4995	0.5005	0.5015	0.5025	0.5035	0.5045	0.5055	0.5065	0.5075	0.5085	0.5095	0.5105	0.5115	0.5125	0.5135	0.5145	0.5155	0.5165	0.5175	0.5185	0.5195	0.5205	0.5215	0.5225	0.5235	0.5245	0.5255	0.5265	0.5275	0.5285	0.5295	0.5305	0.5315	0.5325	0.5335	0.5345	0.5355	0.5365	0.5375	0.5385	0.5395	0.5405	0.5415	0.5425	0.5435	0.5445	0.5455	0.5465	0.5475	0.5485	0.5495	0.5505	0.5515	0.5525	0.5535	0.5545	0.5555	0.5565	0.5575	0.5585	0.5595	0.5605	0.5615	0.5625	0.5635	0.5645	0.5655	0.5665	0.5675	0.5685	0.5695	0.5705	0.5715	0.5725	0.5735	0.5745	0.5755	0.5765	0.5775	0.5785	0.5795	0.5805	0.5815	0.5825	0.5835	0.5845	0.5855	0.5865	0.5875	0.5885	0.5895	0.5905	0.5915	0.5925	0.5935	0.5945	0.5955	0.5965	0.5975	0.5985	0.5995	0.6005	0.6015	0.6025	0.6035	0.6045	0.6055	0.6065	0.6075	0.6085	0.6095	0.6105	0.6115	0.6125	0.6135	0.6145	0.6155	0.6165	0.6175	0.6185	0.6195	0.6205	0.6215	0.6225	0.6235	0.6245	0.6255	0.6265	0.6275	0.6285	0.6295	0.631	0.633	0.635	0.637	0.639	0.641	0.643	0.645	0.647	0.649	0.651	0.653	0.655	0.657	0.659	0.661	0.663	0.665	0.667	0.669	0.671	0.673	0.675	0.677	0.679	0.681	0.683	0.685	0.687	0.689	0.691	0.693	0.695	0.697	0.699	0.701	0.703	0.705	0.707	0.709	0.711	0.713	0.715	0.717	0.719	0.721	0.723	0.725	0.727	0.729	0.731	0.733	0.735	0.737	0.739	0.741	0.743	0.745	0.747	0.749	0.751	0.753	0.755	0.757	0.759	0.761	0.763	0.765	0.767	0.769	0.771	0.773	0.775	0.777	0.779	0.781	0.783	0.785	0.787	0.789	0.791	0.793	0.795	0.797	0.799	0.801	0.803	0.805	0.807	0.809	0.811	0.813	0.815	0.817	0.819	0.821	0.823	0.825	0.826	0.828	0.83	0.832	0.834	0.836	0.838	0.84	0.842	0.844	0.846	0.848	0.85	0.852	0.854	0.856	0.858	0.86	0.862	0.864	0.866	0.868	0.87	0.872	0.874	0.876	0.878	0.88	0.882	0.884	0.886	0.888	0.89	0.892	0.894	0.896	0.898	0.9	0.902	0.904	0.906	0.908	0.91	0.912	0.914	0.916	0.918	0.92	0.922	0.924	0.926	0.928	0.93	0.932	0.934	0.936	0.938	0.94	0.942	0.944	0.946	0.948	0.95	0.952	0.954	0.956	0.958	0.96	0.962	0.964	0.966	0.968	0.97	0.972	0.974	0.976	0.978	0.98	0.982	0.984	0.986	0.988	0.99	0.992	0.994	0.996	0.998	1	1.002	1.004	1.006	1.008	1.01	1.012	1.014	1.016	1.018	1.02	1.022	1.024	1.026	1.028	1.03	1.032	1.034	1.036	1.038	1.04	1.042	1.044	1.046	1.048	1.05	1.052	1.054	1.056	1.058	1.06	1.062	1.064	1.066	1.068	1.07	1.072	1.074	1.076	1.078	1.08	1.082	1.084	1.086	1.088	1.09	1.092	1.094	1.096	1.098	1.1	1.102	1.104	1.106	1.108	1.11	1.112	1.114	1.116	1.118	1.12	1.122	1.124	1.126	1.128	1.13	1.132	1.134	1.136	1.138	1.14	1.142	1.144	1.146	1.148	1.15	1.152	1.154	1.156	1.158	1.16	1.162	1.164	1.166	1.168	1.17	1.172	1.174	1.176	1.178	1.18	1.182	1.184	1.186	1.188	1.19	1.192	1.194	1.196	1.198	1.2	1.202	1.204	1.206	1.208	1.21	1.212	1.214	1.216	1.218	1.22	1.222	1.224	1.226	1.228	1.23	1.232	1.234	1.236	1.238	1.24	1.242	1.244	1.246	1.248	1.25	1.252	1.254	1.256	1.258	1.26	1.262	1.264	1.266	1.268	1.27	1.272	1.274	1.276	1.278	1.28	1.282	1.284	1.286	1.288	1.29	1.292	1.294	1.296	1.298	1.3	1.302	1.304	1.306	1.308	1.31	1.312	1.314	1.316	1.318	1.32	1.322	1.324	1.326	1.328	1.33	1.332	1.334	1.336	1.338	1.34	1.342	1.344	1.346	1.348	1.35	1.352	1.354	1.356	1.358	1.36	1.362	1.364	1.366	1.368	1.37	1.372	1.374	1.376	1.378	1.38	1.382	1.384	1.386	1.388	1.39	1.392	1.394	1.396	1.398	1.4	1.402	1.404	1.406	1.408	1.41	1.412	1.414	1.416	1.418	1.42	1.422	1.424	1.426	1.428	1.43	1.432	1.434	1.436	1.438	1.44	1.442	1.444	1.446	1.448	1.45	1.452	1.454	1.456	1.458	1.46	1.462	1.464	1.466	1.468	1.47	1.472	1.474	1.476	1.478	1.48	1.482	1.484	1.486	1.488	1.49	1.492	1.494	1.496	1.498	1.5	1.502	1.504	1.506	1.508	1.51	1.512	1.514	1.516	1.518	1.52	1.522	1.524	1.526	1.528	1.53	1.532	1.534	1.536	1.538	1.54	1.542	1.544	1.546	1.548	1.55	1.552	1.554	1.556	1.558	1.56	1.562	1.564	1.566	1.568	1.57	1.572	1.574	1.576	1.578	1.58	1.582	1.584	1.586	1.588	1.59	1.592	1.594	1.596	1.598	1.6	1.602	1.604	1.606	1.608	1.61	1.612	1.614	1.616	1.618	1.62	1.622	1.624	1.626	1.628	1.63	1.632	1.634	1.636	1.638	1.64	1.642	1.644	1.646	1.648	1.65	1.652	1.654	1.656	1.658	1.66	1.662	1.664	1.666	1.668	1.67	1.672	1.674	1.676	1.678	1.68	1.682	1.684	1.686	1.688	1.69	1.692	1.694	1.696	1.698	1.7	1.702	1.704	1.706	1.708	1.71	1.712	1.714	1.716	1.718	1.72	1.722	1.724	1.726	1.728	1.73	1.732	1.734	1.736	1.738	1.74	1.742	1.744	1.746	1.748	1.75	1.752	1.754	1.756	1.758	1.76	1.762	1.764	1.766	1.768	1.77	1.772	1.774	1.776	1.778	1.78	1.782	1.784	1.786	1.788	1.79	1.792	1.794	1.796	1.798	1.8	1.802	1.804	1.806	1.808	1.81	1.812	1.814	1.816	1.818	1.82	1.822	1.824	1.826	1.828	1.83	1.832	1.834	1.836	1.838	1.84	1.842	1.844	1.846	1.848	1.85	1.852	1.854	1.856	1.858	1.86	1.862	1.864	1.866	1.868	1.87	1.872	1.874	1.876	1.878	1.88	1.882	1.884	1.886	1.888	1.89	1.892	1.894	1.896	1.898	1.9	1.902	1.904	1.906	1.908	1.91	1.912	1.914	1.916	1.918	1.92	1.922	1.924	1.926	1.928	1.93	1.932	1.934	1.936	1.938	1.94	1.942	1.944	1.946	1.948	1.95	1.952	1.954	1.956	1.958	1.96	1.962	1.964	1.966	1.968	1.97	1.972	1.974	1.976	1.978	1.98	1.982	1.984	1.986	1.988	1.99	1.992	1.994	1.996	1.998	2	2.002	2.004	2.006	2.008	2.01	2.012	2.014	2.016	2.018	2.02	2.022	2.024	2.026	2.028	2.03	2.032	2.034	2.036	2.038	2.04	2.042	2.044	2.046	2.048	2.05	2.052	2.054	2.056	2.058	2.06	2.062	2.064	2.066	2.068	2.07	2.072	2.074	2.076	2.078	2.08	2.082	2.084	2.086	2.088	2.09	2.092	2.094	2.096	2.098	2.1	2.102	2.104	2.106	2.108	2.11	2.112	2.114	2.116	2.118	2.12	2.122	2.124	2.126	2.128	2.13	2.132	2.134	2.136	2.138	2.14	2.142	2.144	2.146	2.148	2.15	2.152	2.154	2.156	2.158	2.16	2.162	2.164	2.166	2.168	2.17	2.172	2.174	2.176	2.178	2.18	2.182	2.184	2.186	2.188	2.19	2.192	2.194	2.196	2.198	2.2	2.202	2.204	2.206	2.208	2.21	2.212	2.214	2.216	2.218	2.22	2.222	2.224	2.226	2.228	2.23	2.232	2.234	2.236	2.238	2.24	2.242	2.244	2.246	2.248	2.25	2.252	2.254	2.256	2.258	2.26	2.262	2.264	2.266	2.268	2.27	2.272	2.274	2.276	2.278	2.28	2.282	2.284	2.286	2.288	2.29	2.292	2.294	2.296	2.298	2.3	2.302	2.304	2.306	2.308	2.31	2.312	2.314	2.316	2.318	2.32	2.322	2.324	2.326	2.328	2.33	2.332	2.334	2.336	2.338	2.34	2.342	2.344	2.346	2.348	2.35	2.352	2.354	2.356	2.358	2.36	2.362	2.364	2.366	2.368	2.37	2.372	2.374	2.376	2.378	2.38	2.382	2.384	2.386	2.388	2.39	2.392	2.394	2.396	2.398	2.4	2.402	2.404	2.406	2.408	2.41	2.412	2.414	2.416	2.418	2.42	2.422	2.424	2.426	2.428	2.43	2.432	2.434	2.436	2.438	2.44	2.442	2.444	2.446	2.448	2.45	2.452	2.454	2.456	2.458	2.46	2.462	2.464	2.466	2.468	2.47	2.472	2.474	2.476	2.478	2.48	2.482	2.484	2.486	2.488	2.49	2.492	2.494	2.496	2.498	2.5	2.52	2.54	2.56	2.58	2.6	2.62	2.64	2.66	2.68	2.7	2.72	2.74	2.76	2.78	2.8	2.82	2.84	2.86	2.88	2.9	2.92	2.94	2.96	2.98	3	3.02	3.04	3.06	3.08	3.1	3.12	3.14	3.16	3.18	3.2	3.22	3.24	3.26	3.28	3.3	3.32	3.34	3.36	3.38	3.4	3.42	3.44	3.46	3.48	3.5	3.52	3.54	3.56	3.58	3.6	3.62	3.64	3.66	3.68	3.7	3.72	3.74	3.76	3.78	3.8	3.82	3.84	3.86	3.88	3.9	3.92	3.94	3.96	3.98	4	4.02	4.04	4.06	4.08	4.1	4.12	4.14	4.16	4.18	4.2	4.22	4.24	4.26	4.28	4.3	4.32	4.34	4.36	4.38	4.4	4.42	4.44	4.46	4.48	4.5	4.52	4.54	4.56	4.58	4.6	4.62	4.64	4.66	4.68	4.7	4.72	4.74	4.76	4.78	4.8	4.82	4.84	4.86	4.88	4.9	4.92	4.94	4.96	4.98	5	5.05	5.1	5.15	5.2	5.25	5.3	5.35	5.4	5.45	5.5	5.55	5.6	5.65	5.7	5.75	5.8	5.85	5.9	5.95	6	6.05	6.1	6.15	6.2	6.25	6.3	6.35	6.4	6.45	6.5	6.55	6.6	6.65	6.7	6.75	6.8	6.85	6.9	6.95	7	7.05	7.1	7.15	7.2	7.25	7.3	7.35	7.4	7.45	7.5	7.55	7.6	7.65	7.7	7.75	7.8	7.85	7.9	7.95	8	8.05	8.1	8.15	8.2	8.25	8.3	8.35	8.4	8.45	8.5	8.55	8.6	8.65	8.7	8.75	8.8	8.85	8.9	8.95	9	9.05	9.1	9.15	9.2	9.25	9.3	9.35	9.4	9.45	9.5	9.55	9.6	9.65	9.7	9.75	9.8	9.85	9.9	9.95	10	11	12	13	14	15	16	17	18	19	20	25	30	35	40	50	60	80	100	120	150	200	250	300	400	1000		];
    am0_phi_per_m2_per_micron = [	6.19E-02	0.5614	4.901	1.184	4.77E-02	3.43E-02	2.88E-02	3.52E-02	2.13E-02	1.73E-02	3.99E-02	0.1206	3.98E-02	4.13E-02	0.168	4.57E-02	3.80E-02	3.09E-02	2.92E-02	3.97E-02	7.56E-02	6.08E-02	4.21E-02	4.68E-02	5.11E-02	5.09E-02	5.54E-02	7.09E-02	8.49E-02	8.20E-02	7.96E-02	8.70E-02	9.27E-02	0.1163	0.1299	0.2059	0.2144	0.1847	0.1717	0.1675	0.1754	0.1934	0.2228	0.2519	0.2841	0.2973	0.4302	0.3989	0.3875	0.4556	0.5877	0.6616	0.688	0.7252	0.7645	0.9067	1.079	1.22	1.403	1.538	1.576	1.831	2.233	2.243	2.244	2.066	2.311	2.7	3.009	3.291	3.569	3.764	4.165	4.113	3.808	5.21	5.427	6.008	6.191	6.187	6.664	7.326	8.023	8.261	9.217	10.25	10.54	11.08	12.65	15.05	21.38	27.92	33.54	31.3	33.15	40.03	36.15	32.27	35.29	44.37	46.92	47.33	39.58	49.65	63.01	58.97	52.29	39.4	39.92	51.95	47.71	52.12	50.97	53.26	44.74	38.97	51.42	48.59	48.44	41.96	44.12	39.56	51.48	70.6	66.53	60.97	49.39	50.4	55.5	45.65	56.38	60.1	46.01	41.55	51.55	59.57	79.3	101.8	125.4	125.1	104	85.51	89.8	103.6	165.8	249.7	252.7	249.4	250.8	243.8	238.9	267.3	224.4	197.4	196.5	132.6	175.1	242.8	233.8	159.3	85.55	94.63	208.3	294.1	313.5	235.3	163.1	322.7	336.3	322.2	472.7	601.3	580.8	521.9	535.5	508.8	553.2	509.6	507.3	465.5	484	420	455.5	489	620.6	602.5	594.8	555.7	615	611.4	496.5	622.4	729.2	655.9	699.9	662.9	633	633.2	773.9	664.9	710.5	805.1	699.5	688.6	661.3	760.8	875.8	979.5	952.7	917.6	1061	1016	965.7	954.9	921.6	958.9	943.4	809.5	841.8	921.5	958.1	1007	923.8	993	950.6	795.7	939.2	926.4	901.7	897.2	889.8	1050	979.5	907.9	1033	1111	1045	912.3	796	693.6	991.1	970.8	878.1	997.8	996.9	1013	1152	1233	1180	1101	1226	1139	1175	1054	920.2	900.4	1062	1085	1282	1327	1066	1202	1082	791.3	684.1	959.7	1008	1007	1004	984.3	1174	1247	1342	1019	582.3	1026	1314	854.5	928.8	1522	1663	1682	1746	1759	1684	1674	1667	1589	1628	1735	1715	1532	1817	1789	1756	1737	1734	1842	1665	1684	1701	1757	1797	1582	1711	1767	1695	1698	1569	1587	1475	1135	1686	1646	1731	1670	1723	1929	1806	1567	1825	1713	1931	1980	1909	1973	1821	1891	2077	1973	2027	2144	2109	1941	1970	1979	2034	2077	2100	1971	2009	2040	2055	2104	2040	1976	2042	1921	2015	1994	1990	1877	2018	2041	1991	2051	2016	1956	2075	2009	2076	2035	2090	2023	2019	1969	1830	1625	1830	1914	1960	2007	1896	1896	1888	2058	1926	2017	2018	1866	1970	1857	1812	1894	1934	1869	1993	1961	1906	1919	1916	1947	1997	1867	1861	1874	1900	1669	1726	1654	1828	1831	1906	1823	1894	1958	1930	1674	1828	1897	1918	1952	1963	1770	1923	1858	1990	1871	1882	1904	1832	1769	1881	1825	1879	1879	1901	1879	1833	1863	1895	1862	1871	1846	1882	1898	1897	1821	1846	1787	1808	1843	1824	1850	1861	1854	1798	1829	1887	1810	1860	1769	1823	1892	1876	1867	1830	1846	1857	1783	1828	1838	1853	1873	1857	1860	1783	1830	1848	1750	1612	1813	1787	1808	1796	1773	1782	1805	1780	1757	1774	1746	1751	1719	1787	1776	1763	1759	1757	1743	1744	1703	1746	1705	1683	1713	1713	1609	1707	1724	1707	1734	1690	1713	1666	1656	1632	1697	1697	1697	1677	1639	1651	1656	1654	1651	1614	1621	1627	1603	1558	1606	1599	1532	1384	1549	1571	1555	1560	1535	1546	1516	1521	1510	1508	1498	1492	1479	1455	1467	1461	1448	1448	1436	1416	1425	1386	1388	1415	1400	1384	1385	1373	1366	1354	1328	1331	1348	1350	1346	1319	1326	1318	1309	1307	1278	1258	1286	1279	1283	1270	1262	1259	1255	1248	1240	1237	1241	1221	1185	1203	1204	1208	1188	1196	1187	1187	1176	1180	1177	1174	1158	1143	1134	1152	1135	1142	1129	1115	1120	1095	1114	1115	1107	1104	1063	1080	1073	1075	1080	1081	1063	1051	1041	1052	1044	1040	1036	1024	1028	1023	966	996.1	878	975.5	1005	996.9	994.9	999.3	886.2	939.5	974.7	983.3	971.3	964	974.9	955.4	951.1	957.9	938.3	944.3	953	939.4	933.2	938.7	933.9	915.8	891.6	928.5	917.6	902.5	891.6	896.7	907.1	900.4	895.1	890.8	863	858.5	861.2	876.9	867.7	865.1	864.1	854.7	858	843.8	825	832.4	837.5	840.7	836.9	831.7	808	808.2	818.8	815.1	808.9	801.3	794.7	796.9	795.9	793.6	781.5	782.5	777.9	774.6	776.4	769.8	766.1	761.5	754.1	756.7	755.6	752.5	751	747.9	746.9	726.1	713.6	733.5	731.3	726.2	721	713.9	710.7	704.1	702.1	705.4	702.7	698.9	693.7	690.5	681.7	684	677.2	676.1	674.6	671.4	660	664.4	662.2	658.6	654.9	655.7	645.1	641.5	643.8	645.9	639.5	631.7	624.1	632.6	627.6	628	627.2	624.7	609.9	618	620.8	610.3	619.9	615.9	584.9	598.3	596.1	604.2	593.2	597.4	594.5	591.6	590.6	584.3	584.4	583.1	581.5	574.1	579.6	576.9	565.5	570	565.3	567.8	563.8	565.8	556.9	553	553.1	551.4	554.8	552.5	548.9	545.8	547.9	545.5	543.5	532	532.5	533.2	530.3	531.2	527.6	531.5	527.3	518.4	519	523.9	515.9	510.3	518.7	507.5	508.5	516.1	514.5	508.4	494.3	500.3	506.8	494.8	503.9	489	488.2	493.3	494.2	493	489.7	487.5	485.4	484.6	481.7	477.1	479.2	475	472.9	471.9	470.3	465.3	464.2	461.9	463.5	463.3	462.4	457.1	457.4	455.1	453.3	453	449.7	447.8	446.7	441.7	445.3	445.2	443.1	445.1	444	435.6	401.4	425.9	432.8	431.4	425.5	425.4	422.3	422.4	418.4	418.6	413.9	411.1	413.6	412.3	410.6	403.3	402.2	397.9	401.7	401.6	398.6	398.1	394.9	390.8	387.8	386.3	389.2	386.6	383.2	379	380.5	379.8	377.2	376.6	372.4	374.2	372.2	367.5	368.8	367.3	367.7	365.7	365.7	362.8	359.9	362.1	361.1	356.1	358	357.9	354.5	354.7	353.2	353	350.6	351.3	348.8	348.7	349.2	342.7	343.9	342.8	343.1	342.7	341.8	334.8	337.7	338.5	338.6	335.7	331.5	331.1	328.1	328.5	325.7	330	328.4	328.5	328.3	318.8	318.6	319.7	321.6	321.6	318.7	315.4	314.3	313.1	316.7	315.6	312.1	310.5	310.8	311.4	310.2	307.3	303.4	304.8	304.4	306.8	304.4	303.9	303.3	285.5	301.5	301.8	303.3	297.2	299.4	301.1	292.4	279.9	284.8	291.9	294.7	291.3	288.3	288.2	288.4	286.6	282.4	283.5	284.6	284.6	276.5	282.3	278.4	280.6	277.3	273	275.3	277.8	277.2	271.1	271.3	273.1	267.6	267.1	268.9	268.3	269.7	266.9	265.4	263.3	264.5	267.3	261	253.6	254.7	265	259	259.1	259.9	249	240.5	252.6	258.3	250.6	254.5	251.2	248.9	249.7	247.7	249.1	240	243	244.9	237.4	242.3	236.9	238.3	241.6	240.2	241.8	239.3	238.7	235.9	235.7	227.4	226.2	226.6	227.8	229.4	229.2	227.2	226.8	226.2	226	225.2	224.5	224.6	222.7	221.2	219.3	222.5	217.3	219.3	216.1	216.8	208	205.4	212.9	213.1	212	210.5	212.3	211.2	210	208.9	206.3	204.7	205.2	205	201.7	201.3	198.2	203.7	202.2	201	199.3	197.5	195.4	198.2	197.1	198.4	193.6	187.4	182.7	186.3	190.5	190.2	190.7	186.7	187.2	185.8	185	185.6	184.9	184.3	183.1	179.3	180.7	181.7	180.2	179.1	179.4	179.2	176.3	174.7	175.6	174.7	173.5	173.9	174.7	173.3	172.1	170.9	170.6	170.3	169.9	167.2	168.8	168.8	168.5	168.6	167.5	165.8	160.5	152	159.6	159.8	162.4	162.8	161.1	160.6	159.3	158.5	158.1	156.2	156.2	154	154.1	153.5	151	154.6	153.4	152.5	150.9	152.5	150.3	150.4	150.9	149.4	149.2	150.8	147.3	140.1	129.9	144.1	146.2	147.4	146.4	143.9	145.3	142.4	140.8	139.6	137.3	139	139.7	140.9	138.6	139	137.7	137.8	135.4	137	136	135.3	133.3	135	134.1	134.4	132.2	131.3	130.8	132	132.8	132.1	129.9	129.4	120.3	119.2	127.1	126.1	125.5	128.6	127.6	127.1	126.1	124	122.2	123.1	124	123.9	121.3	120.8	122.4	119.4	119.6	120.5	119.7	117.8	119.5	119.8	118	116.2	117.3	115.9	117	116.1	114.8	114.7	115.4	114.9	114.5	113.8	113.7	113.4	111.6	110.7	111.6	111.5	110.7	108.6	109.8	109.2	108.3	106.4	107.8	107.6	107.6	107.1	106.3	105.9	104.7	104.6	104.6	104	102.8	102.3	100.5	102.5	101.9	100.3	100.4	100.9	100.6	100	98.78	98.64	97.72	98.52	98.35	97.88	95.67	95.93	95.8	96.2	96.06	95.77	95.59	95.74	95.13	93.96	94.52	94.36	93.31	93.11	92.75	92.75	91.89	92.08	92.25	92.09	92.1	91.55	90.12	91.1	90.83	90.64	90.06	89.39	89.79	89.57	89.13	88.78	88.74	88.42	87.81	86.86	84.56	78.49	83	85.57	85.91	85.92	85.32	84.25	84.97	84.25	84.57	84.65	82.77	83.04	83.77	83.49	83.18	82.99	82.65	82.3	82.11	79.66	79.66	80.8	81.05	80.72	79.94	79.7	79.97	79.62	79.26	78.11	78.26	78.31	78.15	78.02	77.58	76.48	76.39	76.42	76.24	76.12	75.2	75.41	75.12	74.02	74.22	74.41	74.21	72.99	73.29	73.15	73.27	72.97	72.77	72.52	72.39	72.42	71.65	70.07	71.25	71.24	71.27	71.1	70.67	69.2	69.08	69.19	69.53	69.55	69.31	69.23	69.01	68.7	68.67	68.26	67.79	67.45	67.68	66.75	65.36	65.59	66.29	66.16	65.84	65.71	65.36	64.96	65.2	65.39	65.09	64.86	64.72	64.53	62.89	62.39	62.82	62.66	63.08	63.05	62.95	62.84	62.63	62.11	62.07	60.66	61.64	61.92	61.72	60.98	58.85	59.08	60.04	60.29	60.08	60.03	59.96	59.89	59.44	59.65	59.45	59.19	59.15	59.02	58.94	57.34	55.99	57.48	57.7	57.67	57.26	57.17	57.12	57.12	57.02	56.41	56.18	55.99	56.39	56.17	56.03	54.98	54.57	54.62	54.32	54.55	53.7	53.92	54.57	54.42	54.35	54.05	53.9	52.85	53.3	53.13	53.43	53.03	51.77	51.4	52.19	51.6	51.69	52.25	51.98	51.75	51.52	51.54	51.55	49.84	48.14	46.72	45.5	44.57	43.05	42.11	40.79	39.68	38.67	37.63	36.63	35.46	34.68	33.85	32.97	32.09	31.19	30.32	29.69	28.9	28.17	27.5	26.82	26.12	25.47	24.65	24.22	23.64	23.06	22.46	21.98	21.44	20.96	20.48	20	19.51	19.07	18.58	18.02	17.68	17.37	16.97	16.59	16.15	15.84	15.54	15.2	14.86	14.56	14.25	13.93	13.62	13.34	13.07	12.81	12.51	12.22	11.93	11.62	11.45	11.08	10.96	10.78	10.57	10.38	10.19	9.983	9.782	9.599	9.427	9.233	9.032	8.857	8.669	8.557	8.385	8.217	8.054	7.894	7.739	7.587	7.439	7.294	7.153	7.015	6.881	6.749	6.621	6.496	6.374	6.254	6.138	6.024	5.913	5.804	5.698	5.594	5.492	5.393	5.296	5.201	5.108	5.018	4.929	4.842	4.757	4.674	4.593	4.514	4.436	4.36	4.285	4.212	4.141	4.071	4.003	3.936	3.87	3.806	3.743	3.681	3.621	3.562	3.504	3.394	3.267	3.146	3.03	2.92	2.815	2.715	2.619	2.527	2.439	2.355	2.275	2.198	2.124	2.054	1.986	1.921	1.859	1.799	1.742	1.687	1.634	1.583	1.534	1.487	1.442	1.399	1.357	1.317	1.278	1.24	1.204	1.17	1.136	1.104	1.073	1.043	1.014	0.9862	0.9592	0.9331	0.908	0.8836	0.8601	0.8374	0.8154	0.7942	0.7736	0.7537	0.7344	0.7158	0.6977	0.6802	0.6633	0.6469	0.631	0.6156	0.6006	0.5862	0.5721	0.5585	0.5453	0.5324	0.52	0.5079	0.4961	0.4847	0.4737	0.4629	0.4525	0.4423	0.4324	0.4228	0.4135	0.4044	0.3956	0.387	0.3787	0.3706	0.3627	0.355	0.3475	0.3402	0.3331	0.3262	0.3195	0.3129	0.3065	0.3003	0.2942	0.2883	0.2825	0.2769	0.2714	0.2661	0.2608	0.2558	0.2508	0.246	0.2412	0.1635	0.1152	8.34E-02	6.19E-02	4.69E-02	3.61E-02	2.83E-02	2.25E-02	1.81E-02	1.47E-02	6.05E-03	2.93E-03	1.58E-03	9.31E-04	3.83E-04	1.85E-04	5.88E-05	2.42E-05	1.17E-05	4.81E-06	1.53E-06	6.28E-07	2.95E-07	1.01E-07	3.38E-09	];
    % am0 is W/m2*micron
    
    phi_per_nm_vec = am0_phi_per_m2_per_micron / 1e3;
    lambda_nm_vec = lambda_micron_vec * 1e3;
    
    phi_per_nm = interp1(lambda_nm_vec, phi_per_nm_vec, lambda_nm, 'pchip', 'extrap');
    
    q = 1.602e-19; % (C) Electron charge
    E_ev = 1240/lambda_nm;
    E_J = E_ev*q;
    
    flux_per_m2 = phi_per_nm * lambda_nm/E_J;
end