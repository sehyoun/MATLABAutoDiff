function x = exp(x)
    x.values = exp(x.values);
    x.derivatives = valXder(x.values, x.derivatives);