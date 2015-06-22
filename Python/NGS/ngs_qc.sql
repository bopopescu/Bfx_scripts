-- count how many variants for each sample by qscore condition
select count(*),subtype,sample,qscore from gilead_ngsqc group by subtype,sample,qscore order by subtype,sample,qscore;
count | subtype | sample | qscore 
-------+---------+--------+--------
   190 | 1a      | 504    | 22n22c
   187 | 1a      | 504    | 23n23c
   188 | 1a      | 504    | 24n24c
   190 | 1a      | 504    | 25n25c
   188 | 1a      | 504    | 26n26c
   190 | 1a      | 504    | 27n27c
   195 | 1a      | 504    | 28n28c
   192 | 1a      | 504    | 29n29c
   185 | 1a      | 504    | 30n30c
   179 | 1a      | 504    | 31n31c
   156 | 1a      | 504    | 32n32c
   151 | 1a      | 506    | 22n22c
   150 | 1a      | 506    | 23n23c
   150 | 1a      | 506    | 24n24c
   150 | 1a      | 506    | 25n25c
   150 | 1a      | 506    | 26n26c
   149 | 1a      | 506    | 27n27c
   148 | 1a      | 506    | 28n28c
   145 | 1a      | 506    | 29n29c
   139 | 1a      | 506    | 30n30c
   134 | 1a      | 506    | 31n31c
   113 | 1a      | 506    | 32n32c
   212 | 1a      | 516    | 22n22c
   211 | 1a      | 516    | 23n23c
   206 | 1a      | 516    | 24n24c
   206 | 1a      | 516    | 25n25c
   208 | 1a      | 516    | 26n26c
   206 | 1a      | 516    | 27n27c
   207 | 1a      | 516    | 28n28c
   208 | 1a      | 516    | 29n29c
   202 | 1a      | 516    | 30n30c
   196 | 1a      | 516    | 31n31c
   177 | 1a      | 516    | 32n32c
   208 | 1a      | 518    | 22n22c
   206 | 1a      | 518    | 23n23c
   207 | 1a      | 518    | 24n24c
   207 | 1a      | 518    | 25n25c
   208 | 1a      | 518    | 26n26c
   207 | 1a      | 518    | 27n27c
   204 | 1a      | 518    | 28n28c
   202 | 1a      | 518    | 29n29c
   198 | 1a      | 518    | 30n30c
   188 | 1a      | 518    | 31n31c
   165 | 1a      | 518    | 32n32c
   171 | 1b      | 509    | 22n22c
   170 | 1b      | 509    | 23n23c
   165 | 1b      | 509    | 24n24c
   163 | 1b      | 509    | 25n25c
   163 | 1b      | 509    | 26n26c
   160 | 1b      | 509    | 27n27c
   156 | 1b      | 509    | 28n28c
   155 | 1b      | 509    | 29n29c
   145 | 1b      | 509    | 30n30c
   139 | 1b      | 509    | 31n31c
   124 | 1b      | 509    | 32n32c
   138 | 1b      | 510    | 22n22c
   135 | 1b      | 510    | 23n23c
   134 | 1b      | 510    | 24n24c
   130 | 1b      | 510    | 25n25c
   128 | 1b      | 510    | 26n26c
   128 | 1b      | 510    | 27n27c
   128 | 1b      | 510    | 28n28c
   121 | 1b      | 510    | 29n29c
   115 | 1b      | 510    | 30n30c
   107 | 1b      | 510    | 31n31c
    99 | 1b      | 510    | 32n32c
   176 | 1b      | 515    | 22n22c
   173 | 1b      | 515    | 23n23c
   178 | 1b      | 515    | 24n24c
   178 | 1b      | 515    | 25n25c
   175 | 1b      | 515    | 26n26c
   176 | 1b      | 515    | 27n27c
   176 | 1b      | 515    | 28n28c
   179 | 1b      | 515    | 29n29c
   176 | 1b      | 515    | 30n30c
   170 | 1b      | 515    | 31n31c
   166 | 1b      | 515    | 32n32c
(77 rows)

-- Count variants from PVD prediction as control
select count(*),subtype,sample from gilead_ngsqc_pvd group by subtype,sample order by subtype,sample;
 count | subtype | sample 
-------+---------+--------
   104 | 1a      | 504
    69 | 1a      | 506
    44 | 1a      | 516
    67 | 1a      | 518
    86 | 1b      | 509
    97 | 1b      | 510
    86 | 1b      | 515
(7 rows)

-- Look at average coverage by qscore
select avg(cast(coverage as integer)) as avg_coverge,subtype,sample,qscore from gilead_ngsqc group by subtype,sample,qscore order by subtype,sample,qscore;
       avg_coverge      | subtype | sample | qscore 
-----------------------+---------+--------+--------
    30104.742105263158 | 1a      | 504    | 22n22c
    27090.016042780749 | 1a      | 504    | 23n23c
    23943.430851063830 | 1a      | 504    | 24n24c
    22039.142105263158 | 1a      | 504    | 25n25c
    19836.654255319149 | 1a      | 504    | 26n26c
    17614.842105263158 | 1a      | 504    | 27n27c
    14705.538461538462 | 1a      | 504    | 28n28c
    12228.609375000000 | 1a      | 504    | 29n29c
 9967.2648648648648649 | 1a      | 504    | 30n30c
 7414.7988826815642458 | 1a      | 504    | 31n31c
 5166.6987179487179487 | 1a      | 504    | 32n32c
    32424.026490066225 | 1a      | 506    | 22n22c
    28763.486666666667 | 1a      | 506    | 23n23c
    25555.846666666667 | 1a      | 506    | 24n24c
    23358.266666666667 | 1a      | 506    | 25n25c
    20974.540000000000 | 1a      | 506    | 26n26c
    18403.523489932886 | 1a      | 506    | 27n27c
    15654.743243243243 | 1a      | 506    | 28n28c
    12936.468965517241 | 1a      | 506    | 29n29c
    10322.820143884892 | 1a      | 506    | 30n30c
 7601.8134328358208955 | 1a      | 506    | 31n31c
 5222.1504424778761062 | 1a      | 506    | 32n32c
    35377.372641509434 | 1a      | 516    | 22n22c
    31317.369668246445 | 1a      | 516    | 23n23c
    27354.553398058252 | 1a      | 516    | 24n24c
    24815.558252427184 | 1a      | 516    | 25n25c
    21899.706730769231 | 1a      | 516    | 26n26c
    19069.718446601942 | 1a      | 516    | 27n27c
    15793.526570048309 | 1a      | 516    | 28n28c
    12586.908653846154 | 1a      | 516    | 29n29c
 9811.3019801980198020 | 1a      | 516    | 30n30c
 7075.7755102040816327 | 1a      | 516    | 31n31c
 4570.8587570621468927 | 1a      | 516    | 32n32c
    30483.677884615385 | 1a      | 518    | 22n22c
    27214.694174757282 | 1a      | 518    | 23n23c
    23958.623188405797 | 1a      | 518    | 24n24c
    21744.371980676329 | 1a      | 518    | 25n25c
    19262.788461538462 | 1a      | 518    | 26n26c
    16865.202898550725 | 1a      | 518    | 27n27c
    14113.044117647059 | 1a      | 518    | 28n28c
    11542.193069306931 | 1a      | 518    | 29n29c
 8986.9393939393939394 | 1a      | 518    | 30n30c
 6684.7500000000000000 | 1a      | 518    | 31n31c
 4438.0000000000000000 | 1a      | 518    | 32n32c
    30089.421052631579 | 1b      | 509    | 22n22c
    26994.805882352941 | 1b      | 509    | 23n23c
    24095.818181818182 | 1b      | 509    | 24n24c
    22264.288343558282 | 1b      | 509    | 25n25c
    20170.638036809816 | 1b      | 509    | 26n26c
    18124.343750000000 | 1b      | 509    | 27n27c
    15675.512820512821 | 1b      | 509    | 28n28c
    13318.600000000000 | 1b      | 509    | 29n29c
    10882.365517241379 | 1b      | 509    | 30n30c
 8163.8417266187050360 | 1b      | 509    | 31n31c
 5629.3709677419354839 | 1b      | 509    | 32n32c
    34697.442028985507 | 1b      | 510    | 22n22c
    30732.103703703704 | 1b      | 510    | 23n23c
    27175.216417910448 | 1b      | 510    | 24n24c
    25536.376923076923 | 1b      | 510    | 25n25c
    23359.335937500000 | 1b      | 510    | 26n26c
    21042.125000000000 | 1b      | 510    | 27n27c
    18332.929687500000 | 1b      | 510    | 28n28c
    16142.000000000000 | 1b      | 510    | 29n29c
    13490.052173913043 | 1b      | 510    | 30n30c
    10518.869158878505 | 1b      | 510    | 31n31c
 7275.7070707070707071 | 1b      | 510    | 32n32c
    41018.085227272727 | 1b      | 515    | 22n22c
    37158.739884393064 | 1b      | 515    | 23n23c
    32430.382022471910 | 1b      | 515    | 24n24c
    30195.337078651685 | 1b      | 515    | 25n25c
    27692.720000000000 | 1b      | 515    | 26n26c
    25158.926136363636 | 1b      | 515    | 27n27c
    21903.022727272727 | 1b      | 515    | 28n28c
    18603.189944134078 | 1b      | 515    | 29n29c
    15589.420454545455 | 1b      | 515    | 30n30c
    12132.641176470588 | 1b      | 515    | 31n31c
 8363.3132530120481928 | 1b      | 515    | 32n32c

-- Average coverage for PVD
select avg(cast(coverage as integer)) as avg_coverge,subtype,sample from gilead_ngsqc_pvd group by subtype,sample order by subtype,sample;
    avg_coverge     | subtype | sample 
--------------------+---------+--------
 44098.076923076923 | 1a      | 504
 45743.159420289855 | 1a      | 506
 53380.818181818182 | 1a      | 516
 45456.850746268657 | 1a      | 518
 42156.895348837209 | 1b      | 509
 49464.546391752577 | 1b      | 510
 56330.976744186047 | 1b      | 515
(7 rows)

-- Build summary table by position
select distinct
num.sample,
num.con_pos,
n22.var_type as n22c22,
n23.var_type as n23c23,
n24.var_type as n24c24,
n25.var_type as n25c25,
n26.var_type as n26c26,
n27.var_type as n27c27,
n28.var_type as n28c28,
n29.var_type as n29c20,
n30.var_type as n30c30,
n31.var_type as n31c31,
n32.var_type as n32c32

FROM 
(select distinct a.con_pos, a.sample from gilead_ngsqc a where a.sample = '509') num 
  left outer join (select distinct b.con_pos, b.var_type, b.sample from gilead_ngsqc b where b.qscore = '22n22c' ) n22 using (con_pos,sample)
  left outer join (select distinct c.con_pos, c.var_type, c.sample from gilead_ngsqc c where c.qscore = '23n23c' ) n23 using (con_pos,sample)
  left outer join (select distinct d.con_pos, d.var_type, d.sample from gilead_ngsqc d where d.qscore = '24n24c' ) n24 using (con_pos,sample)
  left outer join (select distinct e.con_pos, e.var_type, e.sample from gilead_ngsqc e where e.qscore = '25n25c' ) n25 using (con_pos,sample)
  left outer join (select distinct f.con_pos, f.var_type, f.sample from gilead_ngsqc f where f.qscore = '26n26c' ) n26 using (con_pos,sample)
  left outer join (select distinct g.con_pos, g.var_type, g.sample from gilead_ngsqc g where g.qscore = '27n27c' ) n27 using (con_pos,sample)
  left outer join (select distinct h.con_pos, h.var_type, h.sample from gilead_ngsqc h where h.qscore = '28n28c' ) n28 using (con_pos,sample)
  left outer join (select distinct i.con_pos, i.var_type, i.sample from gilead_ngsqc i where i.qscore = '29n29c' ) n29 using (con_pos,sample)
  left outer join (select distinct j.con_pos, j.var_type, j.sample from gilead_ngsqc j where j.qscore = '30n30c' ) n30 using (con_pos,sample)
  left outer join (select distinct k.con_pos, k.var_type, k.sample from gilead_ngsqc k where k.qscore = '31n31c' ) n31 using (con_pos,sample)
  left outer join (select distinct l.con_pos, l.var_type, l.sample from gilead_ngsqc l where l.qscore = '32n32c' ) n32 using (con_pos,sample)

ORDER BY num.sample,num.con_pos