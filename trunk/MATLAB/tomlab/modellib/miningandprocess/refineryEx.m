% function Result = refineryEx(PriLev)
%
% Creates a TOMLAB LP problem for refinery production
%
% REFINERY
%
% A refinery produces butane, petrol, diesel oil, and heating oil
% from two crudes. Four types of operations are necessary to obtain
% these products: separation, conversion, upgrading, and blending.
%
% The separation phase consists of distilling the raw product into,
% among others, butane, naphtha, gasoil, and a residue. The residue
% subsequently undergoes a conversion phase (catalytic cracking) to
% obtain lighter products. The different products that come out of 
% the distillation are purified (desulfurization or sweetening) or 
% upgraded by a reforming operation that augments their octane value.
% Finally, to obtain the products that will be sold, the refinery
% blends several of the intermediate products in order to fulfill the
% prescribed characteristics of the commercial products. The
% following drawing gives a simplified overview on the production
% processes in this refinery. 
%
% After the distillation, crude 1 gives 3% butane, 15% naphtha, 40%
% gasoil, and 15% residue. Crude 2 results in 5% butane, 20% naphtha,
% 35% gasoil, and 10% residue. The reforming of the naphtha gives 15%
% butane and 85% of reformate (reformed naphtha). The catalytic
% cracking of the residue results in 40% of cracked naphtha and 35%
% of cracked gasoil (note that these percentages do not add up to
% 100% because the process also produces 15% of gas, 5% coke and
% another type of residue that are not taken into consideration in
% our example). The petrol is produced with three ingredients:
% reformed naphtha (reformate), butane, and cracked naphtha. The
% diesel oil is obtained by blending sweetened gasoil, cracked
% gasoil, and cracked naphtha. The heating oil may contain gasoil and
% cracked naphtha without any restrictions on their proportions.
%
%          Crudes
%           |
%           V
%          Distillation -----------------------+
%           |                                  |   
%   +-------+--------+--------------+          |   
%   |                |              |          |   
%  Gasoil           Residue        Naphta      |   
%   |                |              |          |   
%   V                V              V          V   
%  Desulf.       Cat. crack        Reform --> Butane 
%   |             |      |          |          |   
%  Sw.Gasoil  C.Naphta  C.Gasoil Reformate     |   
%   |          |         |          |          V
%   |          |         |          |  To Petrol Blending and Butane
%   |          |         |          V
%   |          |         |        To Petrol Blending
%   |          |         V
%   |          |       To Heating Oil and Diesel Blending
%   |          V
%   |      To Petrol, Heating Oil and Diesel Blending
%   V
%  To Heating Oil and Diesel Blending
%
% Simplified representation of a refinery
%
%
% Certain conditions on the quality of the petrol and diesel oil are
% imposed by law. There are three important characteristics for
% petrol: the octane value, vapor pressure and volatility. The octane
% value is a measure of the anti-knock power in the motor. The vapor
% pressure is a measure of the risk of explosion during storage,
% especially with hot weather. The volatility is a measure for how
% easy the motor is started during cold weather. Finally, the
% maximum sulfur content of the diesel oil is imposed by
% antipollution specifications. The following table summarizes the
% required characteristics of the final products and the composition
% of the intermediate ones. Fields are left empty if no particular
% limit applies. We work with the assumption that all these 
% characteristics blend linearly by weight (in reality, this is only
% true for the sulfur contents).
%
% Characteristics of intermediate and final products
%+--------------+------+---------+-------+-------+-------+------+------+
%|Characteristic|Butane|Reformate|Cracked|Cracked|Desulf |Petrol|Diesel|
%|              |      |         |naphtha|gasoil |gasoil |      | oil  |
%+--------------+------+---------+-------+-------+-------+------+------+
%|Octane value  |  120 |   100   | 74    |    -  |   -   |>=94  |  -   |
%|Vapor pressure|   60 |   2.6   |  4.1  |    -  |   -   |<=12.7|  -   |
%|Volatility    |  105 |   3     | 12    |    -  |   -   |>=17  |  -   |
%|Sulfur (in %) |    - |    -    |  0.12 |  0.76 | 0.03  |  -   |<=0.05|
%+--------------+------+---------+-------+-------+-------+------+------+
%
% In the next month the refinery needs to produce 20,000 tonnes of
% butane, 40,000 tonnes of petrol, 30,000 tonnes of diesel oil, and
% 42,000 tonnes of heating oil. 250,000 tonnes of crude 1 and 500,000
% tonnes of crude 2 are available. The monthly capacity of the
% reformer are 30,000 tonnes, for the desulfurization 40,000 tonnes
% and for the cracking 50,000 tonnes. The cost of processing is based
% on the use of fuel and catalysts for the different operations: the
% costs of distillation, reforming, desulfurization, and cracking are
% $ 2.10, $ 4.18, $ 2.04 and $ 0.60 per tonne respectively.
%
% VARIABLES
%
% demand                     Demand of each product
% supply                     Supply of crudes
% capacity                   Capacity of reformer, desulfer and cracker.
% costs                      Cost for each step.
% supplycomp                 Contents of crudes.
% numvbls                    Number of variables
% octane                     Octane in intermediates
% vappres                    Vapor Pressure in intermediates
% volatility                 Volatility of intermediates
% sulfur                     Sulfur in intermediates
% reformer                   Output of components from reformer
% cracker                    Output of components from cracker
% petrolspec_L               Lower bounds (contents in petrol)
% petrolspec_U               Upper bounds (contents in petrol)
% dieselspec_L               Lower bounds (contents in diesel)
% dieselspec_U               Upper bounds (contents in diesel)
% 
% RESULTS
%
% x_k is a vector containing:
% 4 Final productS                butane, petrol, diesel, heating oil      
% 3 Intermediates to Petrol       petbutane, reformerate, petcrknaphtha 
% 3 Intermediates to Diesel       dslgasoil, dslcrknaphtha, dslcrkgasoil
% 3 Intermediates to Heating Oil  hogasoil, hocrknaphtha, hocrkgasoil   
% 4 Intermediates from Distilling distbutane, naphtha, residue, gasoil  
% 2 Intermediates from Reforming  refbutane, reformerate                
% 2 Intermediates from Cracking   crknaphtha, crkgasoil                 
% 2 Crudes                        crude1, crude2                        
%
% REFERENCES
%
% Applications of optimization... Gueret, Prins, Seveaux
% http://web.univ-ubs.fr/lester/~sevaux/pl/index.html
%
% INPUT PARAMETERS
% PriLev       Print Level
%
% OUTPUT PARAMETERS
% Result       Result structure.

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2006 by Tomlab Optimization Inc., $Release: 5.1.0$
% Written Oct 10, 2005.   Last modified Jan 4, 2006.

function Result = refineryEx(PriLev)

if nargin < 1
   PriLev = 1;
end

demand          = [20;40;30;42]*1000;
supply          = [250;500]*1000;
capacity        = [30;40;50]*1000;
costs           = [2.1;2.1;4.18;.6;2.04];
supplycomp      = [.03 .15 .15 .40; .05 .2 .1 .35];
numvbls         = 4+3+3+3+4+2+2+2;

octane          = [120; 100; 74];
vappres         = [60; 2.6; 4.1];
volatility      = [105; 3; 12];
sulfur          = [.03; .12; .76];

reformer        = [.15 .85]';
cracker         = [.40 .35]';

petrolspec_L    = [94; 17];
petrolspec_U    = [12.7];
dieselspec_L    = [-inf];
dieselspec_U    = [.05];

% PRODUCTS
% butane, petrol, diesel, heating        (FINAL PRODUCTS, 4) 4
% petbutane, reformerate, petcrknaphtha  (IPETROL, 3) 7
% dslgasoil, dslcrknaphtha, dslcrkgasoil (IDIESEL, 3) 10
% hogasoil, hocrknaphtha, hocrkgasoil    (IHO, 3) 13
% distbutane, naphtha, residue, gasoil   (IDIST, 4) 17
% refbutane, reformerate                 (IREF, 2) 19
% crknaphtha, crkgasoil                  (ICRACK, 2) 21
% crude1, crude2                         (CRUDES, 2) 23

Prob = refinery(demand, supply, capacity, costs, supplycomp, numvbls, octane,...
   vappres, volatility, sulfur, reformer, cracker, petrolspec_L, petrolspec_U, ...
   dieselspec_L, dieselspec_U);

Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   disp('4 Final products')
   disp(['   butane        - ' num2str(Result.x_k(1))])
   disp(['   petrol        - ' num2str(Result.x_k(2))])
   disp(['   diesel        - ' num2str(Result.x_k(3))])
   disp(['   heating oil   - ' num2str(Result.x_k(4))])
   disp('3 Intermediates to Petrol')
   disp(['   petbutane     - ' num2str(Result.x_k(5))])
   disp(['   reformerate   - ' num2str(Result.x_k(6))])
   disp(['   petcrknaphtha - ' num2str(Result.x_k(7))])
   disp('3 Intermediates to Diesel')
   disp(['   dslgasoil     - ' num2str(Result.x_k(8))])
   disp(['   dslcrknaphtha - ' num2str(Result.x_k(9))])
   disp(['   dslcrkgasoil  - ' num2str(Result.x_k(10))])
   disp('3 Intermediates to Heating Oil')
   disp(['   hogasoil      - ' num2str(Result.x_k(11))])
   disp(['   hocrknaphtha  - ' num2str(Result.x_k(12))])
   disp(['   hocrkgasoil   - ' num2str(Result.x_k(13))])
   disp('4 Intermediates from Distilling')
   disp(['   distbutane    - ' num2str(Result.x_k(14))])
   disp(['   naphtha       - ' num2str(Result.x_k(15))])
   disp(['   residue       - ' num2str(Result.x_k(16))])
   disp(['   gasoil        - ' num2str(Result.x_k(17))])
   disp('2 Intermediates from Reforming')
   disp(['   refbutane     - ' num2str(Result.x_k(18))])
   disp(['   reformerate   - ' num2str(Result.x_k(19))])          
   disp('2 Intermediates from Cracking')
   disp(['   crknaphtha    - ' num2str(Result.x_k(20))])
   disp(['   crkgasoil     - ' num2str(Result.x_k(21))])         
   disp('2 Crudes')
   disp(['   crude1        - ' num2str(Result.x_k(22))])
   disp(['   crude2        - ' num2str(Result.x_k(23))])               
end

% MODIFICATION LOG
%
% 051007 med   Created
% 060104 med   Updated and corrected
% 060110 per   Added documentation.
% 060125 per   Moved disp to end
