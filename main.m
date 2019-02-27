% Design of Rectangular Footing

clc;
clear all;
format short g;

load value.mat
load input.mat
disp ("Design of Rectangular Footing")
disp("\n")
%Compute the Area of the Footing
Load_on_Column = Load;
printf("Load_on_Column = %d KN \n",Load_on_Column)

Assume_the_Self_Weight_of_Footing = (Self_Weight*100);
printf("Assume_the_Self_Weight_of_Footing = %d %% \n",Assume_the_Self_Weight_of_Footing)

Load_inc = Self_Weight*Load;
Self_Weight_of_Footing = Load_inc;
printf("Self_Weight_of_Footing = %d KN \n",Self_Weight_of_Footing)

Vertical_Load_on_Column = Load+Load_inc;
printf("Vertical_Load_on_Column = %d KN \n",Vertical_Load_on_Column)

disp("\n")
Required_Area_of_Footing = (Vertical_Load_on_Column/Soil_Pressure);
printf("Required_Area_of_Footing = %d m^2 \n",Required_Area_of_Footing)

Provided_Area_of_Footing = round(Vertical_Load_on_Column/Soil_Pressure);
printf("Provided_Area_of_Footing = %d m^2 \n",Provided_Area_of_Footing)

% Side of Rectangular Footing
Length_of_Longer_Side_of_Footing = round(sqrt(1.5*(Provided_Area_of_Footing)));
printf("Length_of_Longer_Side_of_Footing = %d m \n",Length_of_Longer_Side_of_Footing)

Length_of_Shorter_Side_of_Footing = round(Provided_Area_of_Footing/Length_of_Longer_Side_of_Footing);
printf("Length_of_Shorter_Side_of_Footing = %d m \n",Length_of_Shorter_Side_of_Footing)

% Net Upward pressure
Net_Upward_Pressure = (Load)/(Length_of_Longer_Side_of_Footing*Length_of_Shorter_Side_of_Footing);
printf("Net_Upward_Pressure = %d KN/m^2 \n",Net_Upward_Pressure)

% Depth on the basis of moment
disp("\n")
disp("Depth on the basis of Bending Compression")

% Longer Side
M1 = ((Net_Upward_Pressure*Length_of_Shorter_Side_of_Footing)*(Length_of_Longer_Side_of_Footing - Longer_Side_of_Column)*(Length_of_Longer_Side_of_Footing - Longer_Side_of_Column))/8;
Factored_M1 = 1.5 * M1;


% Shorter Side
M2 = M1 = ((Net_Upward_Pressure*Length_of_Longer_Side_of_Footing)*(Length_of_Shorter_Side_of_Footing - Shorter_Side_of_Column)*(Length_of_Shorter_Side_of_Footing - Shorter_Side_of_Column))/8;
Factored_M2 = 1.5 * M2;


xu_max_by_d = (700)/(1100+0.87*Fy);
Ru = (0.36*Fck*xu_max_by_d*(1-0.416*xu_max_by_d));

% Depth on longer side
Depth1 = round((sqrt((Factored_M1*1000000)/(Ru*Length_of_Shorter_Side_of_Footing*1000)))/10)*10;

Depth2 = round((sqrt((Factored_M2*1000000)/(Ru*Length_of_Longer_Side_of_Footing*1000)))/10)*10;

if(Depth1>Depth2)
Depth_on_the_basis_of_Moment = Depth1;
printf("Depth_on_the_basis_of_Moment = %d mm \n",Depth_on_the_basis_of_Moment)
elseif(Depth1<Depth2)
Depth_on_the_basis_of_Moment = Depth2;
printf("Depth_on_the_basis_of_Moment = %d mm \n",Depth_on_the_basis_of_Moment)
endif

disp("\n")
disp("Depth on the Basis of One Way Shear")

Permissible_Shear_Stress = interp2(tables,tables,tables,Fck,pt);
printf("Permissible_Shear_Stress = %d N/mm^2 \n",Permissible_Shear_Stress) 

Depth_on_the_basis_of_One_Way_Shear = round((((0.75*Net_Upward_Pressure*Length_of_Longer_Side_of_Footing-0.75*Net_Upward_Pressure*Longer_Side_of_Column))/(Permissible_Shear_Stress+0.0015*Net_Upward_Pressure))/10)*10;
printf("Depth_on_the_basis_of_One_Way_Shear = %d mm \n",Depth_on_the_basis_of_One_Way_Shear)
disp("\n")
if(Depth_on_the_basis_of_Moment>Depth_on_the_basis_of_One_Way_Shear)
Depth_of_Footing = Depth_on_the_basis_of_Moment;
printf("Depth_of_Footing = %d mm \n",Depth_of_Footing)
elseif(Depth_on_the_basis_of_Moment<Depth_on_the_basis_of_One_Way_Shear)
Depth_of_Footing = Depth_on_the_basis_of_One_Way_Shear;
printf("Depth_of_Footing = %d mm \n",Depth_of_Footing)
endif

% Check for Two Way Shear
disp("\n")
disp("Check for Two Way Shear")
Perimeter = (2*((Longer_Side_of_Column+(Depth_of_Footing/1000))+(Shorter_Side_of_Column+(Depth_of_Footing/1000))))*1000;
printf("Perimeter = %d mm \n",Perimeter)

Side1 = (Longer_Side_of_Column+(Depth_of_Footing/1000));
Side2 = (Shorter_Side_of_Column+(Depth_of_Footing/1000));

Area = Side1*Side2;
printf("Area = %d m^2 \n",Area)

Punching_Shear = 1.5*Net_Upward_Pressure*((Length_of_Longer_Side_of_Footing*Length_of_Shorter_Side_of_Footing)-(Area));
printf("Punching_Shear = %d KN \n",Punching_Shear)

Actual_Shear_Stress = (Punching_Shear*1000)/(Perimeter*Depth_of_Footing);
printf("Actual_Shear_Stress = %d N/mm^2 \n",Actual_Shear_Stress)

Allowable_Shear_Stress = 0.25*sqrt(Fck);
printf("Allowable_Shear_Stress = %d N/mm^2 \n",Allowable_Shear_Stress)

if(Actual_Shear_Stress<Allowable_Shear_Stress)
disp("Hence Safe in Two Way Shear")
elseif(Actual_Shear_Stress>Allowable_Shear_Stress)
disp("Hence not Safe in Two Way Shear")
endif

disp("\n")

Effective_Depth = Depth_of_Footing;
printf("Effective_Depth = %d mm \n",Effective_Depth)

Overall_Depth = Effective_Depth+Clear_Cover;
printf("Overall_Depth = %d mm \n",Overall_Depth)

Effective_Depth_in_one_direction = Overall_Depth-Clear_Cover;
printf("Effective_Depth_in_one_direction = %d mm \n",Effective_Depth_in_one_direction)

Effective_Depth_in_other_direction = Effective_Depth_in_one_direction-dia;
printf("Effective_Depth_in_other_direction = %d mm \n",Effective_Depth_in_other_direction)

disp("\n")


% Design of Reinforcement
disp("Design of Reinforcement for Long Bars")
Ast1 = round((0.5*(Fck/Fy))*(1-sqrt(1-(4.6*Factored_M1*1000000)/(Fck*Length_of_Shorter_Side_of_Footing*1000*Effective_Depth*Effective_Depth)))*Length_of_Shorter_Side_of_Footing*1000*Effective_Depth);
printf("Area_of_Steel = %d mm^2 \n",Ast1)

Area_of_one_bar = ceil((pi/4)*(dia*dia));
printf("Area_of_one_bar = %d mm^2 \n",Area_of_one_bar)

No_of_Bars = ceil(Ast1/Area_of_one_bar)
disp("These are to be distributed uniformly in Length of Shorter Side of Footing")

disp("\n")
disp("Design of Reinforcement for Short Bars")
Ast2 = round((0.5*(Fck/Fy))*(1-sqrt(1-(4.6*Factored_M2*1000000)/(Fck*Length_of_Longer_Side_of_Footing*1000*Effective_Depth*Effective_Depth)))*Length_of_Longer_Side_of_Footing*1000*Effective_Depth);
printf("Area_of_Steel = %d mm^2 \n",Ast2)

Area_of_one_bar = ceil((pi/4)*(dia*dia));
printf("Area_of_one_bar = %d mm^2 \n",Area_of_one_bar)

b = Length_of_Longer_Side_of_Footing/Length_of_Shorter_Side_of_Footing;

Ast2b = (2*Ast2)/(b+1);
printf("Area_of_steel_to_be_provided_in_between_two_dintinct_band = %d mm^2 \n",Ast2b)
No_of_Bars_to_be_provided_in_between_two_dintinct_band = ceil(Ast2b/Area_of_one_bar)

Width_of_each_side_band = 0.5*(Length_of_Longer_Side_of_Footing-Length_of_Shorter_Side_of_Footing);
printf("Width_of_each_side_band = %d m \n",Width_of_each_side_band)

Remaining_Area_in_each_side_band = round(0.5*(Ast2-Ast2b));
printf("Remaining_Area_in_each_side_band = %d mm^2 \n",Remaining_Area_in_each_side_band)

No_of_bars_in_side_band = ceil(Remaining_Area_in_each_side_band/Area_of_one_bar);

if(No_of_bars_in_side_band<3)
No_of_bars_in_side_band = 3
elseif(No_of_bars_in_side_band>3)
No_of_bars_in_side_band = No_of_bars_in_side_band
endif
% Check for Development Lenght
disp("\n")
disp("Check for Development Length")
tbd = interp1 (mat(:,1),mat(:,2),Fck);
Development_Length = round((Fy*dia)/(4*tbd*1.6));
printf("Development_Length = %d mm \n", Development_Length)
Length_Available_of_Bars = round((0.5*(Length_of_Shorter_Side_of_Footing*1000 - Shorter_Side_of_Column*1000) - (Clear_Cover)));
printf("Length_Available_of_Bars = %d mm \n", Length_Available_of_Bars)

if(Development_Length<Length_Available_of_Bars)
disp("Hence Safe in Development Length")
elseif(Development_Length>Length_Available_of_Bars)
disp("Hence not Safe in Development Length")
endif

disp("\n")

% Check for Transfer of Load at the base
disp("Check for Transfer of Load at the base")
A2 = (Longer_Side_of_Column*1000)*(Longer_Side_of_Column*1000);

A1 = ((Longer_Side_of_Column*1000)+2*(2*Overall_Depth))*((Longer_Side_of_Column*1000)+2*(2*Overall_Depth));


square_root_A1_by_A2 = sqrt((A1)/(A2));
if(square_root_A1_by_A2>2)
square_root_A1_by_A2 = 2;
elseif(square_root_A1_by_A2<2)
square_root_A1_by_A2 = square_root_A1_by_A2;
endif

Permissible_Bearing_Stress = 0.45*Fck*square_root_A1_by_A2;
printf("Permissible_Bearing_Stress = %d N/mm^2 \n",Permissible_Bearing_Stress)
Actual_Bearing_Stress = (1.5*Load_on_Column*1000)/(Longer_Side_of_Column*1000*Longer_Side_of_Column*1000);
printf("Actual_Bearing_Stress = %d N/mm^2 \n",Actual_Bearing_Stress)

if(Permissible_Bearing_Stress>Actual_Bearing_Stress)
disp("Hence Safe")
elseif(Permissible_Bearing_Stress<Actual_Bearing_Stress)
disp("Hence not in Safe")
endif
