function inv_A = inv_up( inv_A,      U, inv_C, V,        Add_Remove )


if Add_Remove > 0
     inv_A = inv_A - (inv_A*U) * inv( inv_C + V*inv_A*U ) * (V*inv_A); 
else
     inv_A = inv_A + (inv_A*U) * inv( inv_C - V*inv_A*U ) * (V*inv_A); 
end 