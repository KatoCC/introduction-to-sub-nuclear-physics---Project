//
//  Riportiamo i dati della massa invariante e fittiamo le funzioni di verosimiglianza, applicando le opportune restrizioin ai dati di interesse
//

void massa(){
    // limiti per la massa invarianta della JPsi
    double JPsiMassMin = 2.6;
    double JPsiMassMax = 3.5;
    // numero di bin da plottare nell'istogramma della massa
    int Nbin = 100;
    // creiamo l'istrogramma per la massa
    TH1F *hM = new TH1F("J/#Psi Inv Mass", "J/#Psi Invariant Mass", Nbin, JPsiMassMin, JPsiMassMax);
    // definiamo l'ampiezza dei bin
    float binwidth = (JPsiMassMax - JPsiMassMin) / Nbin;
    // accediamo al file ed inseriamo i dati in un tree
    TFile *file = new TFile("Jpsi_DATA.root");
    TTree *tree = (TTree*)file->Get("data");
    // definiamo le variabili che ci interessa importare dal file
    TLorentzVector *JP = new TLorentzVector;
    int trig, trig1, Jq;    // trigger per vedere se il dato acquisito è un eventuo buono, se ha attivato il trigger, Jq ci indica la carica della particella, che coerentemente deve essere nulla, essento prodotto di muPos e muNeg
    // diamo quindi i giusti indirizzi alle variabili definite in data
    tree->SetBranchAddress("JpsiP", &JP);
    tree->SetBranchAddress("HLT_DoubleMu0_Quarkonium_v1", &trig);
    tree->SetBranchAddress("HLT_DoubleMu0", &trig1);
    tree->SetBranchAddress("JpsiCharge", &Jq);
    // acquisiamo il numero di eventi registrati complessivamente
    int n = tree->GetEntries();
    int n_letti = 0;
    // definiamo quindi il ciclo for sugli eventi acquisiti per selezionare le masse invarianti di interesse
    for(int ev = 0; ev < n; ev++){
        // leggiamo l'evento
        tree->GetEntry(ev);
        // leggiamo la massa invariante
        double JPsiMass = JP->M();
        // definiamo quindi le condizioni di selezione: se le condizioni sono soddisfatte inseriamo l'evento nell'istogramma
        if(JPsiMass > JPsiMassMin && JPsiMass < JPsiMassMax && (trig == 1 || trig1 == 1) && Jq == 0){
            n_letti++;
            hM->Fill(JPsiMass);
        }
    }
    
    TCanvas *c = new TCanvas("c");
    c->cd();
    // definiamo quindi le funzioni di fit per il grafico della massa invariante, in particolare fitteremo il segnale di interesse con una gaussiana, mentre il background con un esponenziale (in particolare, con una coda di un esponenziale)
    TF1 *signal = new TF1("signal", "[0]*exp(-0.5*((x-[1])/[2])**2)", JPsiMassMin, JPsiMassMax);
    TF1 *bkg = new TF1("bkg", "exp([3]+[4]*x)", JPsiMassMin, JPsiMassMax);
    TF1 *tot = new TF1("tot", "[0]*exp(-0.5*((x-[1])/[2])**2) + exp([3]+[4]*x)", JPsiMassMin, JPsiMassMax);
    // forniamo valori iniziali sensati, per la convergenza del fit
    tot->SetParameter(1, 3.096);
    tot->SetParameter(2, 0.05);
    // fittiamo quindi la distribuzione
    hM->Fit("tot", "", "", JPsiMassMin, JPsiMassMax);
    // a partire dal fit, associamo i parametri al segnale e al background
    signal->SetParameters(tot->GetParameter(0),tot->GetParameter(1),tot->GetParameter(2));
    bkg->SetParameters(tot->GetParameter(3), tot->GetParameter(4));
    // calcoliamo l'integrale del segnale
    double Fsig = signal->Integral(JPsiMassMin,JPsiMassMax);
    // calcoliamo l'integrale della funzione totale
    double Ftot = tot->Integral(JPsiMassMin,JPsiMassMax);
    // ?
    cout << "Signal " << Fsig << " : " << tot->GetParameter(0) << "+- " << tot->GetParError(0) <<endl;
    cout << "Backgd " << bkg->Integral(JPsiMassMin,JPsiMassMax)<< endl;
    // vediamo quindi quanti sono gli ingressi complessivi dell'istogramma creato
    double Ntot = hM->GetEntries();
    // calcoliamo il numero di segnali buoni
    double Nsig = Fsig * Ntot / Ftot;
    cout << "Ntot: "<<Ntot<< " FTOT " <<Ftot << " FSIG "<< Fsig << " NSIG " << Nsig<< " bin width " << 1./(binwidth) << endl;
    
    hM->Draw();
    tot->SetLineColor(kGreen+3);
    signal->SetLineColor(kBlue+2);
    bkg->SetLineColor(kRed+1);
    bkg->SetLineStyle(2);
    tot->Draw("same");
    signal->Draw("same");
    bkg->Draw("same");
    // salviamo i plot sia come file root che come immagine png
    c->Print("MassFit.png");
    c->Print("MassFit.root");
}
