function [mcmcparams mcmcrun mcmcruntheta mcmcrunphi]=unpack_mcmcdat(mcmcdat)
%  [mcmcparams mcmcrun mcmcruntheta mcmcrunphi]=unpack_mcmcdat(mcmcdat)
%
% Trivial function to unpack mcmcdat structure
%
% NJB 2012


mcmcparams=mcmcdat.mcmcparams;
mcmcrun=mcmcdat.mcmcrun;
mcmcruntheta=mcmcdat.mcmcruntheta;
mcmcrunphi=mcmcdat.mcmcrunphi;
