FROM rhub/r-minimal:4.5.0

RUN installr -d -a neofetch -t "gfortran git cmake" \
	     codetools caravagnalab/ProCESS@1.1

ARG USER_ID=ProCESS
RUN adduser $USER_ID -D
USER $USER_ID
WORKDIR /home/$USER_ID
ENV HOME=/home/$USER_ID

CMD ["neofetch"]
